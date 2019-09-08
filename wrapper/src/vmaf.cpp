/**
 *
 *  Copyright 2016-2019 Netflix, Inc.
 *
 *     Licensed under the Apache License, Version 2.0 (the "License");
 *     you may not use this file except in compliance with the License.
 *     You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0
 *
 *     Unless required by applicable law or agreed to in writing, software
 *     distributed under the License is distributed on an "AS IS" BASIS,
 *     WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *     See the License for the specific language governing permissions and
 *     limitations under the License.
 *
 */

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <numeric>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <cmath>
#include <iomanip>

#include "libvmaf.h"
#include "vmaf.h"
#include "combo.h"
#include "pugixml/pugixml.hpp"
#include "timer.h"
#include "jsonprint.h"
#include "jsonreader.h"
#include "debug.h"

#define VAL_EQUAL_STR(V,S) (Stringize((V)).compare((S))==0)
#define VAL_IS_LIST(V) ((V).tag=='n') /* check ocval.cc */
#define VAL_IS_NONE(V) ((V).tag=='Z') /* check ocval.cc */
#define VAL_IS_DICT(V) ((V).tag=='t') /* check ocval.cc */

//#define MIN(A,B) ((A)<(B)?(A):(B))

#if !defined(min)
template <class T> static inline T min(T x,T y) { return (x<y)?x:y; }
#endif

static const int BOOTSTRAP_MODEL_NAME_PRECISION = 4;

inline double _round_to_digit(double val, int digit);
std::string _get_file_name(const std::string& s);

inline double _round(double val)
{
    if( val < 0 ) return ceil(val - 0.5);
    return floor(val + 0.5);
}

inline double _round_to_digit(double val, int digit)
{
    size_t m = pow(10.0, digit);
    return _round(val * m) / m;
}

std::string _get_file_name(const std::string& s)
{
    size_t i = s.find_last_of("/\\", s.length());
    if (i != std::string::npos) {
        return(s.substr(i + 1, s.length() - i));
    }
    return("");
}

std::string to_zero_lead(const int value, const unsigned precision)
{
     std::ostringstream oss;
     oss << std::setw(precision) << std::setfill('0') << value;
     return oss.str();
}

std::string pooling_enum_to_string(enum VmafPoolingMethod pool_method) {
    switch (pool_method) {
        case VmafPoolingMethod::VMAF_POOL_MIN: return "min";
        case VmafPoolingMethod::VMAF_POOL_MEAN: return "mean";
        case VmafPoolingMethod::VMAF_POOL_HARMONIC_MEAN: return "harmonic_mean";
        default: return "";
    }
}

void SvmDelete::operator()(void *svm)
{
    svm_free_and_destroy_model((svm_model **)&svm);
}

void LibsvmNusvrTrainTestModel::_assert_model_type(Val model_type) {
    if (!VAL_EQUAL_STR(model_type, "'LIBSVMNUSVR'")) {
        printf("Expect model type LIBSVMNUSVR, "
                "but got %s\n", Stringize(model_type).c_str());
        throw VmafException("Incompatible model_type");
    }
}

void LibsvmNusvrTrainTestModel::_read_and_assert_model(const char *model_path, Val& feature_names,
     Val& norm_type, Val& slopes, Val& intercepts, Val& score_clip, Val& score_transform)
{
    /*  TrainTestModel._assert_trained() in Python:
     *  assert 'model_type' in self.model_dict # need this to recover class
        assert 'feature_names' in self.model_dict
        assert 'norm_type' in self.model_dict
        assert 'model' in self.model_dict

        norm_type = self.model_dict['norm_type']
        assert norm_type == 'none' or norm_type == 'linear_rescale'

        if norm_type == 'linear_rescale':
            assert 'slopes' in self.model_dict
            assert 'intercepts' in self.model_dict
     */
    Val model, model_type;
    char errmsg[1024];
    try
    {
        LoadValFromFile(model_path, model, SERIALIZE_P0);

        model_type = model["model_dict"]["model_type"];
        feature_names = model["model_dict"]["feature_names"];
        norm_type = model["model_dict"]["norm_type"];

        slopes = model["model_dict"]["slopes"];
        intercepts = model["model_dict"]["intercepts"];
        score_clip = model["model_dict"]["score_clip"];
        score_transform = model["model_dict"]["score_transform"];
    }
    catch (std::runtime_error& e)
    {
        printf("Input model at %s cannot be read successfully.\n", model_path);
        sprintf(errmsg, "Error loading model (.pkl): %s", e.what());
        throw VmafException(errmsg);
    }

    _assert_model_type(model_type);

    if (!VAL_IS_LIST(feature_names))
    {
        printf("feature_names in model must be a list.\n");
        throw VmafException("Incompatible feature_names");
    }

    if (!(VAL_EQUAL_STR(norm_type, "'none'") ||
          VAL_EQUAL_STR(norm_type, "'linear_rescale'")
          )
        )
    {
        printf("norm_type in model must be either 'none' or 'linear_rescale'.\n");
        throw VmafException("Incompatible norm_type");
    }

    if ( VAL_EQUAL_STR(norm_type, "'linear_rescale'") &&
         (!VAL_IS_LIST(slopes) || !VAL_IS_LIST(intercepts))
        )
    {
        printf("if norm_type in model is 'linear_rescale', "
                "both slopes and intercepts must be a list.\n");
        throw VmafException("Incompatible slopes or intercepts");
    }

    if (!(VAL_IS_NONE(score_clip) || VAL_IS_LIST(score_clip)))
    {
        printf("score_clip in model must be either None or list.\n");
        throw VmafException("Incompatible score_clip");
    }

    if (!(VAL_IS_NONE(score_transform) || VAL_IS_DICT(score_transform)))
    {
        printf("score_transform in model must be either None or dictionary (table).\n");
        throw VmafException("Incompatible score_transform");
    }

}

void LibsvmNusvrTrainTestModel::load_model()
{
    dbg_printf("Read input model (pkl) at %s ...\n", model_path);
    _read_and_assert_model(model_path, feature_names, norm_type, slopes, intercepts, score_clip, score_transform);

    /* follow the convention that if model_path is a/b.c, the libsvm_model_path is always a/b.c.model */
    std::string libsvm_model_path_ = std::string(model_path) + std::string(".model");
    const char *libsvm_model_path = libsvm_model_path_.c_str();
    dbg_printf("Read input model (libsvm) at %s ...\n", libsvm_model_path);
    svm_model_ptr = _read_and_assert_svm_model(libsvm_model_path);
}

std::unique_ptr<svm_model, SvmDelete> LibsvmNusvrTrainTestModel::_read_and_assert_svm_model(const char* libsvm_model_path)
{
    std::unique_ptr<svm_model, SvmDelete> svm_model_ptr { svm_load_model(libsvm_model_path) };
    if (!svm_model_ptr) {
        printf("Error loading SVM model.\n");
        throw VmafException("Error loading SVM model");
    }
    return svm_model_ptr;
}

void LibsvmNusvrTrainTestModel::populate_and_normalize_nodes_at_frm(size_t i_frm,
        svm_node*& nodes, StatVector& adm2,
        StatVector& adm_scale0, StatVector& adm_scale1,
        StatVector& adm_scale2, StatVector& adm_scale3, StatVector& motion,
        StatVector& vif_scale0, StatVector& vif_scale1, StatVector& vif_scale2,
        StatVector& vif_scale3, StatVector& vif, StatVector& motion2) {
    if (VAL_EQUAL_STR(norm_type, "'linear_rescale'")) {
        for (size_t j = 0; j < feature_names.length(); j++) {
            nodes[j].index = j + 1;
            if (strcmp(Stringize(feature_names[j]).c_str(),
                    "'VMAF_feature_adm2_score'") == 0)
                nodes[j].value = double(slopes[j + 1]) * adm2.at(i_frm)
                        + double(intercepts[j + 1]);
            else if (strcmp(Stringize(feature_names[j]).c_str(),
                    "'VMAF_feature_adm_scale0_score'") == 0)
                nodes[j].value = double(slopes[j + 1])
                        * adm_scale0.at(i_frm)
                        + double(intercepts[j + 1]);
            else if (strcmp(Stringize(feature_names[j]).c_str(),
                    "'VMAF_feature_adm_scale1_score'") == 0)
                nodes[j].value = double(slopes[j + 1])
                        * adm_scale1.at(i_frm)
                        + double(intercepts[j + 1]);
            else if (strcmp(Stringize(feature_names[j]).c_str(),
                    "'VMAF_feature_adm_scale2_score'") == 0)
                nodes[j].value = double(slopes[j + 1])
                        * adm_scale2.at(i_frm)
                        + double(intercepts[j + 1]);
            else if (strcmp(Stringize(feature_names[j]).c_str(),
                    "'VMAF_feature_adm_scale3_score'") == 0)
                nodes[j].value = double(slopes[j + 1])
                        * adm_scale3.at(i_frm)
                        + double(intercepts[j + 1]);
            else if (strcmp(Stringize(feature_names[j]).c_str(),
                    "'VMAF_feature_motion_score'") == 0)
                nodes[j].value = double(slopes[j + 1]) * motion.at(i_frm)
                        + double(intercepts[j + 1]);
            else if (strcmp(Stringize(feature_names[j]).c_str(),
                    "'VMAF_feature_vif_scale0_score'") == 0)
                nodes[j].value = double(slopes[j + 1])
                        * vif_scale0.at(i_frm)
                        + double(intercepts[j + 1]);
            else if (strcmp(Stringize(feature_names[j]).c_str(),
                    "'VMAF_feature_vif_scale1_score'") == 0)
                nodes[j].value = double(slopes[j + 1])
                        * vif_scale1.at(i_frm)
                        + double(intercepts[j + 1]);
            else if (strcmp(Stringize(feature_names[j]).c_str(),
                    "'VMAF_feature_vif_scale2_score'") == 0)
                nodes[j].value = double(slopes[j + 1])
                        * vif_scale2.at(i_frm)
                        + double(intercepts[j + 1]);
            else if (strcmp(Stringize(feature_names[j]).c_str(),
                    "'VMAF_feature_vif_scale3_score'") == 0)
                nodes[j].value = double(slopes[j + 1])
                        * vif_scale3.at(i_frm)
                        + double(intercepts[j + 1]);
            else if (strcmp(Stringize(feature_names[j]).c_str(),
                    "'VMAF_feature_vif_score'") == 0)
                nodes[j].value = double(slopes[j + 1]) * vif.at(i_frm)
                        + double(intercepts[j + 1]);
            else if (strcmp(Stringize(feature_names[j]).c_str(),
                    "'VMAF_feature_motion2_score'") == 0)
                nodes[j].value = double(slopes[j + 1]) * motion2.at(i_frm)
                        + double(intercepts[j + 1]);
            else {
                printf("Unknown feature name: %s.\n",
                        Stringize(feature_names[j]).c_str());
                throw VmafException("Unknown feature name");
            }
        }
    } else {
        for (size_t j = 0; j < feature_names.length(); j++) {
            nodes[j].index = j + 1;
            if (strcmp(Stringize(feature_names[j]).c_str(),
                    "'VMAF_feature_adm2_score'") == 0)
                nodes[j].value = adm2.at(i_frm);
            else if (strcmp(Stringize(feature_names[j]).c_str(),
                    "'VMAF_feature_adm_scale0_score'") == 0)
                nodes[j].value = adm_scale0.at(i_frm);
            else if (strcmp(Stringize(feature_names[j]).c_str(),
                    "'VMAF_feature_adm_scale1_score'") == 0)
                nodes[j].value = adm_scale1.at(i_frm);
            else if (strcmp(Stringize(feature_names[j]).c_str(),
                    "'VMAF_feature_adm_scale2_score'") == 0)
                nodes[j].value = adm_scale2.at(i_frm);
            else if (strcmp(Stringize(feature_names[j]).c_str(),
                    "'VMAF_feature_adm_scale3_score'") == 0)
                nodes[j].value = adm_scale3.at(i_frm);
            else if (strcmp(Stringize(feature_names[j]).c_str(),
                    "'VMAF_feature_motion_score'") == 0)
                nodes[j].value = motion.at(i_frm);
            else if (strcmp(Stringize(feature_names[j]).c_str(),
                    "'VMAF_feature_vif_scale0_score'") == 0)
                nodes[j].value = vif_scale0.at(i_frm);
            else if (strcmp(Stringize(feature_names[j]).c_str(),
                    "'VMAF_feature_vif_scale1_score'") == 0)
                nodes[j].value = vif_scale1.at(i_frm);
            else if (strcmp(Stringize(feature_names[j]).c_str(),
                    "'VMAF_feature_vif_scale2_score'") == 0)
                nodes[j].value = vif_scale2.at(i_frm);
            else if (strcmp(Stringize(feature_names[j]).c_str(),
                    "'VMAF_feature_vif_scale3_score'") == 0)
                nodes[j].value = vif_scale3.at(i_frm);
            else if (strcmp(Stringize(feature_names[j]).c_str(),
                    "'VMAF_feature_vif_score'") == 0)
                nodes[j].value = vif.at(i_frm);
            else if (strcmp(Stringize(feature_names[j]).c_str(),
                    "'VMAF_feature_motion2_score'") == 0)
                nodes[j].value = motion2.at(i_frm);
            else {
                printf("Unknown feature name: %s.\n",
                        Stringize(feature_names[j]).c_str());
                throw VmafException("Unknown feature name");
            }
        }
    }
}

VmafPredictionStruct LibsvmNusvrTrainTestModel::predict(svm_node* nodes) {

    double prediction = svm_predict(svm_model_ptr.get(), nodes);
    VmafPredictionStruct predictionStruct;

    /* denormalize score */
    _denormalize_prediction(prediction);

    predictionStruct.vmafPrediction[VmafPredictionReturnType::SCORE] = prediction;

    return predictionStruct;
}

void LibsvmNusvrTrainTestModel::_denormalize_prediction(double& prediction) {
    if (VAL_EQUAL_STR(norm_type, "'linear_rescale'")) {
        /* denormalize */
        prediction = (prediction - double(intercepts[0]))
                / double(slopes[0]);
    } else {
        ;
    }
}

std::string BootstrapLibsvmNusvrTrainTestModel::_get_model_i_filename(const char* model_path, int i_model)
{
    if (i_model == 0) {
        return std::string(model_path);
    }
    else {
        std::stringstream ss;
        ss << '.' << std::setw(4) << std::setfill('0') << i_model;
        std::string s = ss.str();
        std::string model_path_i = std::string(model_path) + s;
        return model_path_i;
    }
}

void BootstrapLibsvmNusvrTrainTestModel::_assert_model_type(Val model_type) {
    if (!VAL_EQUAL_STR(model_type, "'RESIDUEBOOTSTRAP_LIBSVMNUSVR'") &&
            !VAL_EQUAL_STR(model_type, "'BOOTSTRAP_LIBSVMNUSVR'")) {
        printf("Expect model type BOOTSTRAP_LIBSVMNUSVR or "
                "RESIDUEBOOTSTRAP_LIBSVMNUSVR, but got %s\n", Stringize(model_type).c_str());
        throw VmafException("Incompatible model_type");
    }
}

void BootstrapLibsvmNusvrTrainTestModel::_read_and_assert_model(const char *model_path, Val& feature_names, Val& norm_type, Val& slopes,
        Val& intercepts, Val& score_clip, Val& score_transform, int& numModels)
{
    LibsvmNusvrTrainTestModel::_read_and_assert_model(model_path, feature_names, norm_type, slopes, intercepts, score_clip, score_transform);

    Val model, model_type, numModelsVal;
    char errmsg[1024];
    try
    {
        LoadValFromFile(model_path, model, SERIALIZE_P0);
        numModelsVal = model["param_dict"]["num_models"];
    }
    catch (std::runtime_error& e)
    {
        printf("num_models at %s cannot be read successfully.\n", model_path);
        sprintf(errmsg, "Error loading model (.pkl): %s", e.what());
        throw VmafException(errmsg);
    }
    if (VAL_IS_NONE(numModelsVal))
    {
        printf("num_models cannot be none.\n");
        throw VmafException("num_models cannot be none.");
    }
    numModels = numModelsVal;
}

void BootstrapLibsvmNusvrTrainTestModel::load_model()
{

    int numModels;

    std::string model_path_0 = _get_model_i_filename(this->model_path, 0);

    dbg_printf("Read input model (pkl) at %s ...\n", model_path_0.c_str());
    _read_and_assert_model(model_path_0.c_str(), feature_names, norm_type, slopes, intercepts, score_clip, score_transform, numModels);
    dbg_printf("Number of models: %d\n", numModels);

    for (size_t i=0; i<numModels; i++)
    {
        std::string model_path_i = _get_model_i_filename(this->model_path, i);
        /* follow the convention that if model_path is a/b.c, the libsvm_model_path is always a/b.c.model */
        std::string libsvm_model_path_i = model_path_i + std::string(".model");
        dbg_printf("Read input model (libsvm) at %s ...\n", libsvm_model_path_i.c_str());

        if (i == 0)
        {
            svm_model_ptr = _read_and_assert_svm_model(libsvm_model_path_i.c_str());
        }
        else
        {
            bootstrap_svm_model_ptrs.push_back(_read_and_assert_svm_model(libsvm_model_path_i.c_str()));
        }
    }

}

VmafPredictionStruct BootstrapLibsvmNusvrTrainTestModel::predict(svm_node* nodes) {

    VmafPredictionStruct predictionStruct = LibsvmNusvrTrainTestModel::predict(nodes);

    StatVector bootstrapPredictions;
    double bootstrapPrediction;

    for (size_t i=0; i<bootstrap_svm_model_ptrs.size(); i++) {
        bootstrapPrediction = svm_predict(bootstrap_svm_model_ptrs.at(i).get(), nodes);
        _denormalize_prediction(bootstrapPrediction);
        bootstrapPredictions.append(bootstrapPrediction);
    }

    predictionStruct.vmafPrediction[VmafPredictionReturnType::BAGGING_SCORE] = bootstrapPredictions.mean();
    predictionStruct.vmafPrediction[VmafPredictionReturnType::STDDEV] = bootstrapPredictions.std();
    predictionStruct.vmafPrediction[VmafPredictionReturnType::CI95_LOW] = bootstrapPredictions.percentile(2.5);
    predictionStruct.vmafPrediction[VmafPredictionReturnType::CI95_HIGH] = bootstrapPredictions.percentile(97.5);

    predictionStruct.vmafMultiModelPrediction = bootstrapPredictions.getVector();

    return predictionStruct;
}

void VmafQualityRunner::_transform_value(LibsvmNusvrTrainTestModel& model,
        double& prediction) {
    if (!VAL_IS_NONE(model.score_transform)) {
        double value = 0.0;

        /* quadratic transform */
        if (!VAL_IS_NONE(model.score_transform["p0"])) {
            value += double(model.score_transform["p0"]);
        }
        if (!VAL_IS_NONE(model.score_transform["p1"])) {
            value += double(model.score_transform["p1"]) * prediction;
        }
        if (!VAL_IS_NONE(model.score_transform["p2"])) {
            value += double(model.score_transform["p2"]) * prediction
                    * prediction;
        }

        /* rectification */
        if (!VAL_IS_NONE(model.score_transform["out_lte_in"])
                && VAL_EQUAL_STR(model.score_transform["out_lte_in"], "'true'")) {
            if (value > prediction) {
                value = prediction;
            }
        }
        if (!VAL_IS_NONE(model.score_transform["out_gte_in"])
                && VAL_EQUAL_STR(model.score_transform["out_gte_in"], "'true'")) {
            if (value < prediction) {
                value = prediction;
            }
        }

        prediction = value;
    }
}

void VmafQualityRunner::_clip_value(LibsvmNusvrTrainTestModel& model,
        double& prediction) {
    if (!VAL_IS_NONE(model.score_clip)) {
        if (prediction < double(model.score_clip[0])) {
            prediction = double(model.score_clip[0]);
        } else if (prediction > double(model.score_clip[1])) {
            prediction = double(model.score_clip[1]);
        }
    }
}

void VmafQualityRunner::_postproc_predict(
        VmafPredictionStruct& predictionStruct) {
    ;
}

void VmafQualityRunner::_transform_score(LibsvmNusvrTrainTestModel& model,
        VmafPredictionStruct& predictionStruct) {

    double& prediction = predictionStruct.vmafPrediction[VmafPredictionReturnType::SCORE];
    _transform_value(model, prediction);
}

void VmafQualityRunner::_clip_score(LibsvmNusvrTrainTestModel& model,
        VmafPredictionStruct& predictionStruct) {

    double& prediction = predictionStruct.vmafPrediction[VmafPredictionReturnType::SCORE];
    _clip_value(model, prediction);
}

void VmafQualityRunner::_postproc_transform_clip(
        VmafPredictionStruct& predictionStruct) {
    ;
}

void VmafQualityRunner::_normalize_predict_denormalize_transform_clip(
        LibsvmNusvrTrainTestModel& model, size_t num_frms, Result& result,
        int vmaf_model_setting,
        std::vector<VmafPredictionStruct>& predictionStructs) {

    StatVector adm2 = result.get_scores("adm2");
    StatVector adm_scale0 = result.get_scores("adm_scale0");
    StatVector adm_scale1 = result.get_scores("adm_scale1");
    StatVector adm_scale2 = result.get_scores("adm_scale2");
    StatVector adm_scale3 = result.get_scores("adm_scale3");

    StatVector vif = result.get_scores("vif");
    StatVector vif_scale0 = result.get_scores("vif_scale0");
    StatVector vif_scale1 = result.get_scores("vif_scale1");
    StatVector vif_scale2 = result.get_scores("vif_scale2");
    StatVector vif_scale3 = result.get_scores("vif_scale3");

    StatVector motion = result.get_scores("motion");
    StatVector motion2 = result.get_scores("motion2");

    /* IMPORTANT: always allocate one more spot and put a -1 at the last one's
     * index, so that libsvm will stop looping when seeing the -1 !!!
     * see https://github.com/cjlin1/libsvm */
    svm_node* nodes = (svm_node*) alloca(
            sizeof(svm_node) * (model.feature_names.length() + 1));
    nodes[model.feature_names.length()].index = -1;

    size_t i_subsampled;
    for (size_t i_frm=0; i_frm<num_frms; i_frm++) {
        model.populate_and_normalize_nodes_at_frm(i_frm, nodes, adm2,
                adm_scale0, adm_scale1, adm_scale2, adm_scale3, motion,
                vif_scale0, vif_scale1, vif_scale2, vif_scale3, vif, motion2);

        VmafPredictionStruct predictionStruct = model.predict(nodes);

        _postproc_predict(predictionStruct);

        if (vmaf_model_setting & VMAF_MODEL_SETTING_ENABLE_TRANSFORM) {
            _transform_score(model, predictionStruct);
        }

        if (!(vmaf_model_setting & VMAF_MODEL_SETTING_DISABLE_CLIP)) {
            _clip_score(model, predictionStruct);
        }

        _postproc_transform_clip(predictionStruct);

        dbg_printf("frame: %zu, ", i_frm);
        dbg_printf("adm2: %f, ", adm2.at(i_frm));
        dbg_printf("adm_scale0: %f, ", adm_scale0.at(i_frm));
        dbg_printf("adm_scale1: %f, ", adm_scale1.at(i_frm));
        dbg_printf("adm_scale2: %f, ", adm_scale2.at(i_frm));
        dbg_printf("adm_scale3: %f, ", adm_scale3.at(i_frm));
        dbg_printf("motion: %f, ", motion.at(i_frm));
        dbg_printf("vif_scale0: %f, ", vif_scale0.at(i_frm));
        dbg_printf("vif_scale1: %f, ", vif_scale1.at(i_frm));
        dbg_printf("vif_scale2: %f, ", vif_scale2.at(i_frm));
        dbg_printf("vif_scale3: %f, ", vif_scale3.at(i_frm));
        dbg_printf("vif: %f, ", vif.at(i_frm));
        dbg_printf("motion2: %f, ", motion2.at(i_frm));

        dbg_printf("\n");

        predictionStructs.push_back(predictionStruct);

    }
}

std::unique_ptr<LibsvmNusvrTrainTestModel> VmafQualityRunner::_load_model(const char *model_path)
{
    std::unique_ptr<LibsvmNusvrTrainTestModel> model_ptr = std::unique_ptr<LibsvmNusvrTrainTestModel>(new LibsvmNusvrTrainTestModel(model_path));
    model_ptr->load_model();
    return model_ptr;
}

void VmafQualityRunner::_set_prediction_result(
        std::vector<VmafPredictionStruct> predictionStructs,
        Result& result, std::string model_name) {
    StatVector score;
    for (size_t i = 0; i < predictionStructs.size(); i++) {
        score.append(predictionStructs.at(i).vmafPrediction[VmafPredictionReturnType::SCORE]);
    }
    result.set_scores(model_name, score);
}

void VmafQualityRunner::feature_extract(Result &result,
                        int (*read_vmaf_picture)(VmafPicture *ref_vmaf_pict, VmafPicture *dis_vmaf_pict, float *temp_data, void *user_data),
                        void *user_data, VmafSettings *vmafSettings)
{
    dbg_printf("Initialize storage arrays...\n");

    unsigned int w = vmafSettings->width;
    unsigned int h = vmafSettings->height;
    enum VmafPixelFormat fmt = vmafSettings->pix_fmt;

    char errmsg[1024];
    DArray adm_num_array, adm_den_array, adm_num_scale0_array,
            adm_den_scale0_array, adm_num_scale1_array, adm_den_scale1_array,
            adm_num_scale2_array, adm_den_scale2_array, adm_num_scale3_array,
            adm_den_scale3_array, motion_array, motion2_array,
            vif_num_scale0_array, vif_den_scale0_array, vif_num_scale1_array,
            vif_den_scale1_array, vif_num_scale2_array, vif_den_scale2_array,
            vif_num_scale3_array, vif_den_scale3_array, vif_array,
            psnr_array, psnr_u_array, psnr_v_array,
            ssim_array, ms_ssim_array;
    /* use the following ptrs as flags to turn on/off optional metrics */
    DArray *psnr_array_ptr, *psnr_u_array_ptr, *psnr_v_array_ptr;
    DArray *ssim_array_ptr, *ms_ssim_array_ptr;
    init_array(&adm_num_array, INIT_FRAMES);
    init_array(&adm_den_array, INIT_FRAMES);
    init_array(&adm_num_scale0_array, INIT_FRAMES);
    init_array(&adm_den_scale0_array, INIT_FRAMES);
    init_array(&adm_num_scale1_array, INIT_FRAMES);
    init_array(&adm_den_scale1_array, INIT_FRAMES);
    init_array(&adm_num_scale2_array, INIT_FRAMES);
    init_array(&adm_den_scale2_array, INIT_FRAMES);
    init_array(&adm_num_scale3_array, INIT_FRAMES);
    init_array(&adm_den_scale3_array, INIT_FRAMES);
    init_array(&motion_array, INIT_FRAMES);
    init_array(&motion2_array, INIT_FRAMES);
    init_array(&vif_num_scale0_array, INIT_FRAMES);
    init_array(&vif_den_scale0_array, INIT_FRAMES);
    init_array(&vif_num_scale1_array, INIT_FRAMES);
    init_array(&vif_den_scale1_array, INIT_FRAMES);
    init_array(&vif_num_scale2_array, INIT_FRAMES);
    init_array(&vif_den_scale2_array, INIT_FRAMES);
    init_array(&vif_num_scale3_array, INIT_FRAMES);
    init_array(&vif_den_scale3_array, INIT_FRAMES);
    init_array(&vif_array, INIT_FRAMES);
    init_array(&psnr_array, INIT_FRAMES);
    init_array(&psnr_u_array, INIT_FRAMES);
    init_array(&psnr_v_array, INIT_FRAMES);
    init_array(&ssim_array, INIT_FRAMES);
    init_array(&ms_ssim_array, INIT_FRAMES);

    /* optional output arrays */
    if (vmafSettings->vmaf_feature_mode_setting & VMAF_FEATURE_MODE_SETTING_DO_PSNR) {
        psnr_array_ptr = &psnr_array;
    } else {
        psnr_array_ptr = NULL;
    }
    if (vmafSettings->vmaf_feature_mode_setting & VMAF_FEATURE_MODE_SETTING_DO_SSIM) {
        ssim_array_ptr = &ssim_array;
    } else {
        ssim_array_ptr = NULL;
    }
    if (vmafSettings->vmaf_feature_mode_setting & VMAF_FEATURE_MODE_SETTING_DO_MS_SSIM) {
        ms_ssim_array_ptr = &ms_ssim_array;
    } else {
        ms_ssim_array_ptr = NULL;
    }

    bool use_color = vmafSettings->vmaf_feature_mode_setting & VMAF_FEATURE_MODE_SETTING_DO_COLOR;

    if (vmafSettings->vmaf_feature_mode_setting & VMAF_FEATURE_MODE_SETTING_DO_COLOR) {
        psnr_u_array_ptr = &psnr_u_array;
        psnr_v_array_ptr = &psnr_v_array;
    } else {
        psnr_u_array_ptr = NULL;
        psnr_v_array_ptr = NULL;
    }

    dbg_printf("Extract atom features...\n");
    int ret = combo(read_vmaf_picture, user_data, w, h, fmt, &adm_num_array,
            &adm_den_array, &adm_num_scale0_array, &adm_den_scale0_array,
            &adm_num_scale1_array, &adm_den_scale1_array, &adm_num_scale2_array,
            &adm_den_scale2_array, &adm_num_scale3_array, &adm_den_scale3_array,
            &motion_array, &motion2_array, &vif_num_scale0_array,
            &vif_den_scale0_array, &vif_num_scale1_array, &vif_den_scale1_array,
            &vif_num_scale2_array, &vif_den_scale2_array, &vif_num_scale3_array,
            &vif_den_scale3_array, &vif_array,
            psnr_array_ptr, psnr_u_array_ptr, psnr_v_array_ptr, ssim_array_ptr,
            ms_ssim_array_ptr, errmsg, vmafSettings->vmaf_feature_calculation_setting,
            use_color);
    if (ret) {
        throw VmafException(errmsg);
    }

    size_t num_frms = motion_array.used;
    bool num_frms_is_consistent = (adm_num_array.used == num_frms)
            && (adm_den_array.used == num_frms)
            && (motion2_array.used == num_frms)
            && (adm_num_scale0_array.used == num_frms)
            && (adm_den_scale0_array.used == num_frms)
            && (adm_num_scale1_array.used == num_frms)
            && (adm_den_scale1_array.used == num_frms)
            && (adm_num_scale2_array.used == num_frms)
            && (adm_den_scale2_array.used == num_frms)
            && (adm_num_scale3_array.used == num_frms)
            && (adm_den_scale3_array.used == num_frms)
            && (vif_num_scale0_array.used == num_frms)
            && (vif_den_scale0_array.used == num_frms)
            && (vif_num_scale1_array.used == num_frms)
            && (vif_den_scale1_array.used == num_frms)
            && (vif_num_scale2_array.used == num_frms)
            && (vif_den_scale2_array.used == num_frms)
            && (vif_num_scale3_array.used == num_frms)
            && (vif_den_scale3_array.used == num_frms)
            && (vif_array.used == num_frms);
    if (psnr_array_ptr != NULL) {
        num_frms_is_consistent = num_frms_is_consistent
                && (psnr_array.used == num_frms);
    }
    if (ssim_array_ptr != NULL) {
        num_frms_is_consistent = num_frms_is_consistent
                && (ssim_array.used == num_frms);
    }
    if (ms_ssim_array_ptr != NULL) {
        num_frms_is_consistent = num_frms_is_consistent
                && (ms_ssim_array.used == num_frms);
    }
    if (!num_frms_is_consistent) {
        sprintf(errmsg,
                "Output feature vectors are of inconsistent dimensions: motion (%zu), motion2 (%zu), adm_num (%zu), adm_den (%zu), vif_num_scale0 (%zu), vif_den_scale0 (%zu), vif_num_scale1 (%zu), vif_den_scale1 (%zu), vif_num_scale2 (%zu), vif_den_scale2 (%zu), vif_num_scale3 (%zu), vif_den_scale3 (%zu), vif (%zu), psnr (%zu), ssim (%zu), ms_ssim (%zu), adm_num_scale0 (%zu), adm_den_scale0 (%zu), adm_num_scale1 (%zu), adm_den_scale1 (%zu), adm_num_scale2 (%zu), adm_den_scale2 (%zu), adm_num_scale3 (%zu), adm_den_scale3 (%zu)",
                motion_array.used, motion2_array.used, adm_num_array.used,
                adm_den_array.used, vif_num_scale0_array.used,
                vif_den_scale0_array.used, vif_num_scale1_array.used,
                vif_den_scale1_array.used, vif_num_scale2_array.used,
                vif_den_scale2_array.used, vif_num_scale3_array.used,
                vif_num_scale3_array.used, vif_array.used, psnr_array.used,
                ssim_array.used, ms_ssim_array.used, adm_num_scale0_array.used,
                adm_den_scale0_array.used, adm_num_scale1_array.used,
                adm_den_scale1_array.used, adm_num_scale2_array.used,
                adm_den_scale2_array.used, adm_num_scale3_array.used,
                adm_num_scale3_array.used);
        throw VmafException(errmsg);
    }
    dbg_printf(
            "Generate final features (including derived atom features)...\n");
    double ADM2_CONSTANT = 0.0;
    double ADM_SCALE_CONSTANT = 0.0;
    StatVector adm2, motion, vif_scale0, vif_scale1, vif_scale2, vif_scale3,
            vif, motion2;
    StatVector adm_scale0, adm_scale1, adm_scale2, adm_scale3;
    StatVector psnr, psnr_u, psnr_v, ssim, ms_ssim;

    unsigned int num_frms_subsampled = 0;
    for (size_t i = 0; i < num_frms; i += vmafSettings->vmaf_feature_calculation_setting.n_subsample) {
        adm2.append(
                (get_at(&adm_num_array, i) + ADM2_CONSTANT)
                        / (get_at(&adm_den_array, i) + ADM2_CONSTANT));
        adm_scale0.append(
                (get_at(&adm_num_scale0_array, i) + ADM_SCALE_CONSTANT)
                        / (get_at(&adm_den_scale0_array, i) + ADM_SCALE_CONSTANT));
        adm_scale1.append(
                (get_at(&adm_num_scale1_array, i) + ADM_SCALE_CONSTANT)
                        / (get_at(&adm_den_scale1_array, i) + ADM_SCALE_CONSTANT));
        adm_scale2.append(
                (get_at(&adm_num_scale2_array, i) + ADM_SCALE_CONSTANT)
                        / (get_at(&adm_den_scale2_array, i) + ADM_SCALE_CONSTANT));
        adm_scale3.append(
                (get_at(&adm_num_scale3_array, i) + ADM_SCALE_CONSTANT)
                        / (get_at(&adm_den_scale3_array, i) + ADM_SCALE_CONSTANT));
        motion.append(get_at(&motion_array, i));
        motion2.append(get_at(&motion2_array, i));
        vif_scale0.append(
                get_at(&vif_num_scale0_array, i)
                        / get_at(&vif_den_scale0_array, i));
        vif_scale1.append(
                get_at(&vif_num_scale1_array, i)
                        / get_at(&vif_den_scale1_array, i));
        vif_scale2.append(
                get_at(&vif_num_scale2_array, i)
                        / get_at(&vif_den_scale2_array, i));
        vif_scale3.append(
                get_at(&vif_num_scale3_array, i)
                        / get_at(&vif_den_scale3_array, i));
        vif.append(get_at(&vif_array, i));

        if (psnr_array_ptr != NULL) {
            psnr.append(get_at(&psnr_array, i));
        }
        if (psnr_u_array_ptr != NULL) {
            psnr_u.append(get_at(&psnr_u_array, i));
        }
        if (psnr_v_array_ptr != NULL) {
            psnr_v.append(get_at(&psnr_v_array, i));
        }
        if (ssim_array_ptr != NULL) {
            ssim.append(get_at(&ssim_array, i));
        }
        if (ms_ssim_array_ptr != NULL) {
            ms_ssim.append(get_at(&ms_ssim_array, i));
        }

        num_frms_subsampled += 1;
    }

    result = {};
    result.set_scores("adm2", adm2);
    result.set_scores("adm_scale0", adm_scale0);
    result.set_scores("adm_scale1", adm_scale1);
    result.set_scores("adm_scale2", adm_scale2);
    result.set_scores("adm_scale3", adm_scale3);
    result.set_scores("motion", motion);
    result.set_scores("vif_scale0", vif_scale0);
    result.set_scores("vif_scale1", vif_scale1);
    result.set_scores("vif_scale2", vif_scale2);
    result.set_scores("vif_scale3", vif_scale3);
    result.set_scores("vif", vif);
    result.set_scores("motion2", motion2);

    result.set_num_frms(num_frms_subsampled);

    if (psnr_array_ptr != NULL) {
        result.set_scores("psnr", psnr);
    }
    if (psnr_u_array_ptr != NULL) {
        result.set_scores("psnr_u", psnr_u);
    }
    if (psnr_v_array_ptr != NULL) {
        result.set_scores("psnr_v", psnr_v);
    }
    if (ssim_array_ptr != NULL) {
        result.set_scores("ssim", ssim);
    }
    if (ms_ssim_array_ptr != NULL) {
        result.set_scores("ms_ssim", ms_ssim);
    }

    free_array(&adm_num_array);
    free_array(&adm_den_array);
    free_array(&adm_num_scale0_array);
    free_array(&adm_den_scale0_array);
    free_array(&adm_num_scale1_array);
    free_array(&adm_den_scale1_array);
    free_array(&adm_num_scale2_array);
    free_array(&adm_den_scale2_array);
    free_array(&adm_num_scale3_array);
    free_array(&adm_den_scale3_array);
    free_array(&motion_array);
    free_array(&motion2_array);
    free_array(&vif_num_scale0_array);
    free_array(&vif_den_scale0_array);
    free_array(&vif_num_scale1_array);
    free_array(&vif_den_scale1_array);
    free_array(&vif_num_scale2_array);
    free_array(&vif_den_scale2_array);
    free_array(&vif_num_scale3_array);
    free_array(&vif_den_scale3_array);
    free_array(&vif_array);
    free_array(&psnr_array);
    free_array(&psnr_u_array);
    free_array(&psnr_v_array);
    free_array(&ssim_array);
    free_array(&ms_ssim_array);

}

void VmafQualityRunner::predict(Result &result, VmafModel *vmaf_model_ptr)
{
    dbg_printf("Normalize features, SVM regression, denormalize score, clip...\n");

    int num_frms_subsampled = result.get_num_frms();

    // load model
    std::unique_ptr<LibsvmNusvrTrainTestModel> model_ptr = _load_model(vmaf_model_ptr->path);
    LibsvmNusvrTrainTestModel& model = *model_ptr;

    // verify that all feature names in the loaded model are valid
    for (size_t j = 0; j < model.feature_names.length(); j++) {

        if ((strcmp(Stringize(model.feature_names[j]).c_str(),
                "'VMAF_feature_adm2_score'") != 0)
                && (strcmp(Stringize(model.feature_names[j]).c_str(),
                "'VMAF_feature_adm_scale0_score'") != 0)
                && (strcmp(Stringize(model.feature_names[j]).c_str(),
                "'VMAF_feature_adm_scale1_score'") != 0)
                && (strcmp(Stringize(model.feature_names[j]).c_str(),
                "'VMAF_feature_adm_scale2_score'") != 0)
                && (strcmp(Stringize(model.feature_names[j]).c_str(),
                "'VMAF_feature_adm_scale3_score'") != 0)
                && (strcmp(Stringize(model.feature_names[j]).c_str(),
                "'VMAF_feature_motion_score'") != 0)
                && (strcmp(Stringize(model.feature_names[j]).c_str(),
                "'VMAF_feature_motion2_score'") != 0)
                && (strcmp(Stringize(model.feature_names[j]).c_str(),
                "'VMAF_feature_vif_scale0_score'") != 0)
                && (strcmp(Stringize(model.feature_names[j]).c_str(),
                "'VMAF_feature_vif_scale1_score'") != 0)
                && (strcmp(Stringize(model.feature_names[j]).c_str(),
                "'VMAF_feature_vif_scale2_score'") != 0)
                && (strcmp(Stringize(model.feature_names[j]).c_str(),
                "'VMAF_feature_vif_scale3_score'") != 0)
                && (strcmp(Stringize(model.feature_names[j]).c_str(),
                "'VMAF_feature_vif_score'") != 0)) {
                printf("Unknown feature name: %s.\n", Stringize(model.feature_names[j]).c_str());
                throw VmafException("Unknown feature name");
        }

    }

    std::vector<VmafPredictionStruct> predictionStructs;

    _normalize_predict_denormalize_transform_clip(model, num_frms_subsampled, result,
        vmaf_model_ptr->vmaf_model_setting, predictionStructs);

    _set_prediction_result(predictionStructs, result, vmaf_model_ptr->name);

}

std::unique_ptr<LibsvmNusvrTrainTestModel> BootstrapVmafQualityRunner::_load_model(const char *model_path)
{
    std::unique_ptr<LibsvmNusvrTrainTestModel> model_ptr = std::unique_ptr<BootstrapLibsvmNusvrTrainTestModel>(new BootstrapLibsvmNusvrTrainTestModel(model_path));
    model_ptr->load_model();
    return model_ptr;
}

void BootstrapVmafQualityRunner::_postproc_predict(
        VmafPredictionStruct& predictionStruct) {

    double baggingScore = predictionStruct.vmafPrediction[VmafPredictionReturnType::BAGGING_SCORE];
    double scorePlusDelta = baggingScore + BootstrapVmafQualityRunner::DELTA;
    double scoreMinusDelta = baggingScore - BootstrapVmafQualityRunner::DELTA;

    predictionStruct.vmafPrediction[VmafPredictionReturnType::PLUS_DELTA] = scorePlusDelta;
    predictionStruct.vmafPrediction[VmafPredictionReturnType::MINUS_DELTA] = scoreMinusDelta;
}

void BootstrapVmafQualityRunner::_transform_score(LibsvmNusvrTrainTestModel& model,
        VmafPredictionStruct& predictionStruct) {

    double& score = predictionStruct.vmafPrediction[VmafPredictionReturnType::SCORE];
    double& baggingScore = predictionStruct.vmafPrediction[VmafPredictionReturnType::BAGGING_SCORE];
    double& ci95LowScore = predictionStruct.vmafPrediction[VmafPredictionReturnType::CI95_LOW];
    double& ci95HighScore = predictionStruct.vmafPrediction[VmafPredictionReturnType::CI95_HIGH];
    double& scorePlusDelta = predictionStruct.vmafPrediction[VmafPredictionReturnType::PLUS_DELTA];
    double& scoreMinusDelta = predictionStruct.vmafPrediction[VmafPredictionReturnType::MINUS_DELTA];

    _transform_value(model, score);
    _transform_value(model, baggingScore);
    _transform_value(model, ci95LowScore);
    _transform_value(model, ci95HighScore);
    _transform_value(model, scorePlusDelta);
    _transform_value(model, scoreMinusDelta);

    // transform bootstrap model scores
    size_t num_models = predictionStruct.vmafMultiModelPrediction.size();

    for (size_t model_i = 0; model_i < num_models; model_i++)
    {
        _transform_value(model, predictionStruct.vmafMultiModelPrediction.at(model_i));
    }
}

void BootstrapVmafQualityRunner::_clip_score(LibsvmNusvrTrainTestModel& model,
        VmafPredictionStruct& predictionStruct) {

    double& score = predictionStruct.vmafPrediction[VmafPredictionReturnType::SCORE];
    double& baggingScore = predictionStruct.vmafPrediction[VmafPredictionReturnType::BAGGING_SCORE];
    double& ci95LowScore = predictionStruct.vmafPrediction[VmafPredictionReturnType::CI95_LOW];
    double& ci95HighScore = predictionStruct.vmafPrediction[VmafPredictionReturnType::CI95_HIGH];
    double& scorePlusDelta = predictionStruct.vmafPrediction[VmafPredictionReturnType::PLUS_DELTA];
    double& scoreMinusDelta = predictionStruct.vmafPrediction[VmafPredictionReturnType::MINUS_DELTA];

    _clip_value(model, score);
    _clip_value(model, baggingScore);
    _clip_value(model, ci95LowScore);
    _clip_value(model, ci95HighScore);
    _clip_value(model, scorePlusDelta);
    _clip_value(model, scoreMinusDelta);

    // clip bootstrap model scores
    size_t num_models = predictionStruct.vmafMultiModelPrediction.size();

    for (size_t model_i = 0; model_i < num_models; model_i++)
    {
        _clip_value(model, predictionStruct.vmafMultiModelPrediction.at(model_i));
    }
}

void BootstrapVmafQualityRunner::_postproc_transform_clip(
        VmafPredictionStruct& predictionStruct) {

    double scorePlusDelta = predictionStruct.vmafPrediction[VmafPredictionReturnType::PLUS_DELTA];
    double scoreMinusDelta = predictionStruct.vmafPrediction[VmafPredictionReturnType::MINUS_DELTA];
    double& scoreStdDev = predictionStruct.vmafPrediction[VmafPredictionReturnType::STDDEV];

    double slope = (scorePlusDelta - scoreMinusDelta) / (2.0 * BootstrapVmafQualityRunner::DELTA);
    scoreStdDev *= slope;
}

void BootstrapVmafQualityRunner::_set_prediction_result(
        std::vector<VmafPredictionStruct> predictionStructs,
        Result& result, std::string model_name) {

    VmafQualityRunner::_set_prediction_result(predictionStructs, result, model_name);

    StatVector baggingScore, stdDev, ci95LowScore, ci95HighScore;
    for (size_t i = 0; i < predictionStructs.size(); i++) {
        baggingScore.append(predictionStructs.at(i).vmafPrediction[VmafPredictionReturnType::BAGGING_SCORE]);
        stdDev.append(predictionStructs.at(i).vmafPrediction[VmafPredictionReturnType::STDDEV]);
        ci95LowScore.append(predictionStructs.at(i).vmafPrediction[VmafPredictionReturnType::CI95_LOW]);
        ci95HighScore.append(predictionStructs.at(i).vmafPrediction[VmafPredictionReturnType::CI95_HIGH]);
    }
    result.set_scores(model_name + "_bagging", baggingScore);
    result.set_scores(model_name + "_stddev", stdDev);
    result.set_scores(model_name + "_ci95_low", ci95LowScore);
    result.set_scores(model_name + "_ci95_high", ci95HighScore);

    // num_models is same across frames, so just use first frame length
    size_t num_models = 0; 
    if (predictionStructs.size() > 0) {
        num_models = predictionStructs.at(0).vmafMultiModelPrediction.size();
    }

    std::vector<double> perModelScore;

    // bootstrap model result key example: vmaf_bootstrap_0001 is the first bootstrapped model for vmaf
    for (size_t j = 0; j < num_models; j++) {
        for (size_t i = 0; i < predictionStructs.size(); i++) {
            perModelScore.push_back(predictionStructs.at(i).vmafMultiModelPrediction.at(j));
        }
        result.set_scores(model_name + BOOSTRAP_VMAF_MODEL_KEY + to_zero_lead(j + 1, BOOTSTRAP_MODEL_NAME_PRECISION), perModelScore);
        perModelScore.clear();
    }

}

double RunVmaf(int (*read_vmaf_picture)(VmafPicture *ref_vmaf_pict, VmafPicture *dis_vmaf_pict, float *temp_data, void *user_data),
               void *user_data, VmafSettings *vmafSettings)
{
    printf("Start calculating VMAF score...\n");

    if (vmafSettings->width <= 0)
    {
        throw VmafException("Invalid width value (must be > 0)");
    }
    if (vmafSettings->height <= 0)
    {
        throw VmafException("Invalid height value (must be > 0)");
    }
    if (vmafSettings->vmaf_feature_calculation_setting.n_threads < 0)
    {
        throw VmafException("Invalid n_thread value (must be >= 0)");
    }
    if (vmafSettings->vmaf_feature_calculation_setting.n_subsample <= 0)
    {
        throw VmafException("Invalid n_subsample value (must be > 0)");
    }

    Timer timer;
    timer.start();
    Result result;

    // feature extraction
    VmafQualityRunner::feature_extract(result, read_vmaf_picture, user_data, vmafSettings);

    unsigned int num_models = vmafSettings->num_models;
    for (int i = 0; i < num_models; i ++)
    {
        std::unique_ptr<IVmafQualityRunner> runner_ptr =
            VmafQualityRunnerFactory::createVmafQualityRunner(&(vmafSettings->vmaf_model[i]));
        // predict using i-th model
        runner_ptr->predict(result, &(vmafSettings->vmaf_model[i]));
    }

    timer.stop();

    if (vmafSettings->pool_method != NULL && (vmafSettings->pool_method == VMAF_POOL_MIN))
    {
        result.setScoreAggregateMethod(ScoreAggregateMethod::MINIMUM);
    }
    else if (vmafSettings->pool_method != NULL && (vmafSettings->pool_method == VMAF_POOL_HARMONIC_MEAN))
    {
        result.setScoreAggregateMethod(ScoreAggregateMethod::HARMONIC_MEAN);
    }
    else // mean or default
    {
        result.setScoreAggregateMethod(ScoreAggregateMethod::MEAN);
    }

    unsigned int default_model_ind = vmafSettings->default_model_ind;
    std::vector<std::string> result_keys = result.get_keys();
    std::string default_model_name = vmafSettings->vmaf_model[default_model_ind].name;

    size_t num_frames_subsampled = result.get_scores(default_model_name).size();
    double aggregate_vmaf = result.get_score(default_model_name);

    double exec_fps = (double)num_frames_subsampled * vmafSettings->vmaf_feature_calculation_setting.n_subsample / (double)timer.elapsed();
#if TIME_TEST_ENABLE
	double time_taken = (double)timer.elapsed();
#endif
    printf("Exec FPS: %f\n", exec_fps);

    double aggregate_bagging = 0.0, aggregate_stddev = 0.0, aggregate_ci95_low = 0.0, aggregate_ci95_high = 0.0;
    if (result.has_scores(default_model_name + "_bagging"))
        aggregate_bagging = result.get_score(default_model_name + "_bagging");
    if (result.has_scores(default_model_name + "_stddev"))
        aggregate_stddev = result.get_score(default_model_name + "_stddev");
    if (result.has_scores(default_model_name + "_ci95_low"))
        aggregate_ci95_low = result.get_score(default_model_name + "_ci95_low");
    if (result.has_scores(default_model_name + "_ci95_high"))
        aggregate_ci95_high = result.get_score(default_model_name + "_ci95_high");

    double aggregate_psnr = 0.0, aggregate_ssim = 0.0, aggregate_ms_ssim = 0.0;
    if (result.has_scores("psnr"))
        aggregate_psnr = result.get_score("psnr");
    if (result.has_scores("ssim"))
        aggregate_ssim = result.get_score("ssim");
    if (result.has_scores("ms_ssim"))
        aggregate_ms_ssim = result.get_score("ms_ssim");

    std::string pool_method_str = pooling_enum_to_string(vmafSettings->pool_method);

    printf("VMAF score (%s) = %f\n", pool_method_str.c_str(), aggregate_vmaf);
    if (aggregate_bagging)
        printf("VMAF Bagging score (%s) = %f\n", pool_method_str.c_str(), aggregate_bagging);
    if (aggregate_stddev)
        printf("VMAF StdDev score (%s) = %f\n", pool_method_str.c_str(), aggregate_stddev);
    if (aggregate_ci95_low)
        printf("VMAF CI95_low score (%s) = %f\n", pool_method_str.c_str(), aggregate_ci95_low);
    if (aggregate_ci95_high)
        printf("VMAF CI95_high score (%s) = %f\n", pool_method_str.c_str(), aggregate_ci95_high);

    if (aggregate_psnr)
        printf("PSNR score (%s) = %f\n", pool_method_str.c_str(), aggregate_psnr);
    if (aggregate_ssim)
        printf("SSIM score (%s) = %f\n", pool_method_str.c_str(), aggregate_ssim);
    if (aggregate_ms_ssim)
        printf("MS-SSIM score (%s) = %f\n", pool_method_str.c_str(), aggregate_ms_ssim);

    // print out results for all additional (non-default) models
    for (unsigned int model_ind = 0; model_ind < num_models; model_ind++)
    {
        if (model_ind != default_model_ind) {
            printf("%s score (%s) = %f\n", vmafSettings->vmaf_model[model_ind].name,
                pool_method_str.c_str(), result.get_score(vmafSettings->vmaf_model[model_ind].name));
        }
    }

    int num_bootstrap_models = 0;
    std::string bootstrap_model_list_str = "";

    // determine number of bootstrap models (if any) and construct a comma-separated string of bootstrap vmaf model names
    for (size_t j=0; j<result_keys.size(); j++)
    {
        if (result_keys[j].find(BOOSTRAP_VMAF_MODEL_KEY)!= std::string::npos)
        {
            if (num_bootstrap_models == 0)
            {
                bootstrap_model_list_str += result_keys[j] + ",";
            }
            else if (num_bootstrap_models == 1)
            {
                bootstrap_model_list_str += result_keys[j];
            }
            else
            {
                bootstrap_model_list_str += "," + result_keys[j];
            }
            num_bootstrap_models += 1;
        }
    }

    if (vmafSettings->log_path != NULL && vmafSettings->log_fmt != NULL && (vmafSettings->log_fmt == VMAF_LOG_FMT_JSON))
    {
        /* output to json */

        double value;

        OTab params;
        params["model"] = _get_file_name(std::string(vmafSettings->vmaf_model[default_model_ind].path));
        params["scaledWidth"] = vmafSettings->width;
        params["scaledHeight"] = vmafSettings->height;
        params["subsample"] = vmafSettings->vmaf_feature_calculation_setting.n_subsample;
        params["num_bootstrap_models"] = num_bootstrap_models;
        params["bootstrap_model_list_str"] = bootstrap_model_list_str;

        Arr metrics;
        for (size_t j=0; j<result_keys.size(); j++)
        {
            metrics.append(result_keys[j]);
        }

        Arr frames;
        for (size_t i_subsampled=0; i_subsampled<num_frames_subsampled; i_subsampled++)
        {
            OTab frame;
            frame["frameNum"] = i_subsampled * vmafSettings->vmaf_feature_calculation_setting.n_subsample;
            OTab metrics_scores;
            for (size_t j=0; j<result_keys.size(); j++)
            {
                value = result.get_scores(result_keys[j].c_str()).at(i_subsampled);
                value = _round_to_digit(value, 5);
                metrics_scores[result_keys[j].c_str()] = value;
            }
            frame["metrics"] = metrics_scores;
            frames.append(frame);
        }

        Val top = OTab();
        top["version"] = VMAFOSS_DOC_VERSION;
        top["params"] = params;
        top["metrics"] = metrics;
        top["frames"] = frames;

        top["VMAF score"] = aggregate_vmaf;
        top["ExecFps"] = exec_fps;
        if (aggregate_psnr)
            top["PSNR score"] = aggregate_psnr;
        if (aggregate_ssim)
            top["SSIM score"] = aggregate_ssim;
        if (aggregate_ms_ssim)
            top["MS-SSIM score"] = aggregate_ms_ssim;

        std::ofstream log_file(vmafSettings->log_path);
        JSONPrint(top, log_file, 0, true, 2);
        log_file.close();
    }
	else if (vmafSettings->log_path != NULL && vmafSettings->log_fmt != NULL && (vmafSettings->log_fmt == VMAF_LOG_FMT_CSV))
	{
		/* output to csv */
		FILE *csv = fopen(vmafSettings->log_path, "wt");
		fprintf(csv, "Frame,Width,Height,");
		for (size_t j = 0; j < result_keys.size(); j++)
		{
			fprintf(csv, "%s,", result_keys[j].c_str());
		}
		fprintf(csv, "\n");
		for (size_t i_subsampled = 0; i_subsampled<num_frames_subsampled; i_subsampled++)
		{
			int frameNum = i_subsampled * vmafSettings->vmaf_feature_calculation_setting.n_subsample;
			fprintf(csv, "%d,%d,%d,", frameNum, vmafSettings->width, vmafSettings->height);
			for (size_t j = 0; j<result_keys.size(); j++)
			{
				fprintf(csv, "%4.4f,", (float)result.get_scores(result_keys[j].c_str()).at(i_subsampled));
			}
			fprintf(csv, "\n");
		}
		fclose(csv);
	}
    else if ((vmafSettings->log_path != NULL) || (vmafSettings->log_fmt == VMAF_LOG_FMT_XML))
    {
        /* output to xml */

        pugi::xml_document xml;
        pugi::xml_node xml_root = xml.append_child("VMAF");
        xml_root.append_attribute("version") = VMAFOSS_DOC_VERSION;

        auto params_node = xml_root.append_child("params");
        params_node.append_attribute("model") = _get_file_name(std::string(vmafSettings->vmaf_model[default_model_ind].path)).c_str();
        params_node.append_attribute("scaledWidth") = vmafSettings->width;
        params_node.append_attribute("scaledHeight") = vmafSettings->height;
        params_node.append_attribute("subsample") = vmafSettings->vmaf_feature_calculation_setting.n_subsample;
        params_node.append_attribute("num_bootstrap_models") = num_bootstrap_models;
        params_node.append_attribute("bootstrap_model_list_str") = bootstrap_model_list_str.c_str();

        auto info_node = xml_root.append_child("fyi");
        info_node.append_attribute("numOfFrames") = (int)num_frames_subsampled;
        info_node.append_attribute("aggregateVMAF") = aggregate_vmaf;
        if (aggregate_bagging)
            info_node.append_attribute("aggregateBagging") = aggregate_bagging;
        if (aggregate_stddev)
            info_node.append_attribute("aggregateStdDev") = aggregate_stddev;
        if (aggregate_ci95_low)
            info_node.append_attribute("aggregateCI95_low") = aggregate_ci95_low;
        if (aggregate_ci95_high)
            info_node.append_attribute("aggregateCI95_high") = aggregate_ci95_high;
        if (aggregate_psnr)
            info_node.append_attribute("aggregatePSNR") = aggregate_psnr;
        if (aggregate_ssim)
            info_node.append_attribute("aggregateSSIM") = aggregate_ssim;
        if (aggregate_ms_ssim)
            info_node.append_attribute("aggregateMS_SSIM") = aggregate_ms_ssim;
        if (vmafSettings->pool_method)
            info_node.append_attribute("poolMethod") = pool_method_str.c_str();
        info_node.append_attribute("execFps") = exec_fps;
#if TIME_TEST_ENABLE
		info_node.append_attribute("timeTaken") = time_taken;
#endif

        auto frames_node = xml_root.append_child("frames");
        for (size_t i_subsampled=0; i_subsampled<num_frames_subsampled; i_subsampled++)
        {
            auto node = frames_node.append_child("frame");
            node.append_attribute("frameNum") = (int)i_subsampled * vmafSettings->vmaf_feature_calculation_setting.n_subsample;
            for (size_t j=0; j<result_keys.size(); j++)
            {
                node.append_attribute(result_keys[j].c_str()) = result.get_scores(result_keys[j].c_str()).at(i_subsampled);
            }
        }

        xml.save_file(vmafSettings->log_path);
    }

    return aggregate_vmaf;
}
