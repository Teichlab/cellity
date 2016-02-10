## testing function for assessing cell quality with the SVM approach

context("tests on assessment of cell quality with SVM")

test_that("the assess_cell_quality_SVM() function works on example data", {
    data(param_mES_all)
    data(training_mES_features)
    data(training_mES_labels)
    data(mES1_features)
    data(mES1_labels)
    mES1_features_all <- mES1_features[[1]]
    training_mES_features_all <- training_mES_features[[1]]
    mES1_quality_SVM <- assess_cell_quality_SVM(
        training_mES_features_all, training_mES_labels[,2], param_mES_all, 
        mES1_features_all)
    
    expect_that(mES1_quality_SVM, is_a("data.frame"))
})

