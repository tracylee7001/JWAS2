
    input_parameters = InputParameters(options.chainLength     # number of iterations
                                       options.probFixed       # parameter "pi" the probability SNP effect is zero
                                       options.estimatePi      # "yes" or "no"
                                       options.estimateScale   # "yes" or "no"
                                       options.dfEffectVar     # hyper parameter (degrees of freedom) for locus effect variance
                                       options.nuRes           # hyper parameter (degrees of freedom) for residual variance
                                       options.varGenotypic    # used to derive hyper parameter (scale) for locus effect variance
                                       options.varResidual     # used to derive hyper parameter (scale) for locus effect variance
                                       )

    marker_matrix = MarkerMatrix(X)
    fixed_matrix  = FixedMatrix(C)
    estimated_parameters = EstimatedParameters(input_parameters,marker_matrix,fixed_matrix)
    output= OutputValues(input_parameters,marker_matrix,fixed_matrix)
    data=Data(y,Rinv)


function MCMC(input_parameters::InputParameters,
              marker_matrix::MarkerMatrix,
              fixed_matrix::FixedMatrix,
              estimated_parameters::EstimatedParameters,
              output::OutputValues)

    for i=1:input_parameters.chainlength
        #sample fixed effects
        sampleFixedEffects!(yCorr,fixed_matrix,input_parameters,output)

        # sample residula variance
        output.resVar[i] = sampleVariance(yCorr.*RinvSqrt, marker_matrix.nObs,
                                          input_parameters.nuRes,
                                          EstimatedParameters.scaleRes)

        # sample locus effect variance
        for j=1:nMarkers
            locusEffectVar[j] = sampleVariance(α[j],1,dfEffectVar, scaleVar)
        end
