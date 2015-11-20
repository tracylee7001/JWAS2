include("Tools.jl")
include("Samplers.jl")

abstract BayesMethod
type BayesB <: BayesMethod
type BayesC <: BayesMethod
type BayesA <: BayesB

type InputParameters
  method::AbstractString
  chainLength::Int64     # number of iterations
  probFixed::Float64     # parameter "pi" the probability SNP effect is zero
  estimatePi::bool       # "yes" or "no"
  estimateScale::bool    # "yes" or "no"
  dfEffectVar::Float64   # hyper parameter (degrees of freedom) for locus effect variance
  nuRes::Float64         # hyper parameter (degrees of freedom) for residual variance
  varGenotypic::Float64  # used to derive hyper parameter (scale) for locus effect variance
  varResidual::Float64   # used to derive hyper parameter (scale) for locus effect variance
end

type OutputValues
    meanFxdEff::Array{Float64,1}
    meanMrkEff::Array{Float64,2}
    mdlFrq::Array{Float64,1}
    resVar::Array{Float64,1}
    genVar::Array{Float64,1}
    pi::Array{Float64,1}
    scale::Array{Float64,1}

    function OutputValues(input_parameters,marker_matrix,fixed_matrix)
      nFixedEffects = fixed_matrix.nFixedEffects
      nMarkers      = marker_matrix.nMarkers
      chainLength   = input_parameters.chainLength
      new(zeros(nFixedEffects),
          zeros(nMarkers,1),
          zeros(nMarkers),
          zeros(chainLength),
          zeros(chainLength),
          zeros(chainLength),
          zeros(chainLength))
    end
end

type EstimatedParameters
    vare::Float64
    varEffects::Float64
    scaleVar::Float64   # scale factor for locus effects
    scaleRes::Float64   # scale factor for residual varianc
    β                   # sample of fixed effects
    α                   # sample of partial marker effects unconditional on δ
    u                   # sample of marker effects
    δ                   # inclusion indicator for marker effects
    π                   # probFixed
    locusEffectVar      #locuda-specific variance

    function EstimatedParameters(input_parameters::InputParameters,
                                 marker_matrix::MarkerMatrix
                                 fixed_matrix::FixedMatrix)
        vare        = input_parameters.varResidual
        π           = input_parameters.probFixed
        dfEffectVar = input_parameters.dfEffectVar
        nuRes       = input_parameters.nuRes
        varGenotypic= input_parameters.varGenotypic
        mean2pq     = marker_matrix.mean2pq
        nMarkers    = marker_matrix.nMarkers
        nFixedEffects  = fixed_matrix.nFixedEffects

        varEffects     = varGenotypic/((1-π)*mean2pq)
        locusEffectVar = fill(varEffects,nMarkers)
        scaleVar       = varEffects*(dfEffectVar-2)/dfEffectVar # scale factor for locus effects
        scaleRes       = vare*(nuRes-2)/nuRes                   # scale factor for residual varianc
        β          = zeros(nFixedEffects)  # sample of fixed effects
        α          = zeros(nMarkers)       # sample of partial marker effects unconditional on δ
        δ          = zeros(nMarkers)       # inclusion indicator for marker effects
        u          = zeros(nMarkers)       # sample of marker effects
        new(vare,varEffects,scaleVar,scaleRes,β,α,u,δ,π,locusEffectVar)
    end
end

type MarkerMatrix
    X::Array{Float64,2}
    nObs::Int64         #maybe put them into a seperate type all numbers
    nMarkers::Int64
    centered::Bool
    alleleFreq::Array{Float64,2}
    mean2pq::Float64
    xArray::Array{Array{Float64,1},1}
    XpRinvX::Array{Float64,1}

    function MarkerMatrix(X::Array{Float64,2})
        nObs,nMarkers = size(X)
        markerMeans   = center!(X) #centering
        centered      = true
        p             = markerMeans/2.0
        mean2pq       = (2*p*(1-p)')[1,1]
        xArray        = get_column_ref(X)
        XpRinvX       = getXpRinvX(X) ##maybe outside
        new(X,nObs,nMarkers,centered,p,mean2pq,xArray,XpRinvX)
    end
end

type FixedMatrix
    C::Array{Float64,2}
    nFixedEffects = size(C,2)
    Rinv::Array{Float64,2}
    cArray::Array{Array{Float64,1},1}
    CpRinvC::Array{Float64,1}
    boolRinv::Bool
    fixedNames::Array{AbstractString,1}

    function FixedMatrix(file::DataFrames) ###More
        mean2pq     = (2*p*(1-p)')[1,1]
        xArray = get_column_ref(X)
        XpRinvX = getXpRinvX(X)
        new(X,xArray,XpRinvX,markerMeans,mean2pq,true)
    end
end

type Data
      y::Array{Float64,1}
      yCorr::Array{Float64,1}
      Rinv::Array{Float64}
end

type myJWAS
  method::BayesMethod
  input::InputParameters
  output::OutputValues
  marker_matrix::MarkerMatrix
  fixed_matrix::Array{Float64,2}
end


#################################

    #initial values
    vare       = varResidual
    p          = markerMeans/2.0
    mean2pq    = (2*p*(1-p)')[1,1]
    varEffects = varGenotypic/((1-probFixed)*mean2pq)
    scaleVar   = varEffects*(dfEffectVar-2)/dfEffectVar        # scale factor for locus effects
    scaleRes   = varResidual*(nuRes-2)/nuRes        # scale factor for residual varianc
    yCorr      = copy(y)
    β          = zeros(nFixedEffects)  # sample of fixed effects
    α          = zeros(nMarkers)       # sample of partial marker effects unconditional on δ
    u          = zeros(nMarkers)       # sample of marker effects
    δ          = zeros(nMarkers)       # inclusion indicator for marker effects
    π          = probFixed
    RinvSqrt   = sqrt(Rinv)
    locusEffectVar = fill(varEffects,nMarkers)

    #return values
    meanFxdEff = zeros(nFixedEffects)
    meanMrkEff  = zeros(nMarkers,1)
    mdlFrq     = zeros(nMarkers)
    resVar     = zeros(chainLength)
    genVar     = zeros(chainLength)
    pi         = zeros(chainLength)
    scale      = zeros(chainLength)

#################################

    numIter         =   chainLength
    nObs,nMarkers = size(X)
    nFixedEffects = size(C,2)

    ###START
    X = copy(X) #copy original X matrix, not change original data
    markerMeans = center!(X)
    xArray = get_column_ref(X)
    XpRinvX = getXpRinvX(X, Rinv)

    #initial values
    vare       = varResidual
    p          = markerMeans/2.0
    mean2pq    = (2*p*(1-p)')[1,1]
    varEffects = varGenotypic/((1-probFixed)*mean2pq)
    scaleVar   = varEffects*(dfEffectVar-2)/dfEffectVar        # scale factor for locus effects
    scaleRes   = varResidual*(nuRes-2)/nuRes        # scale factor for residual varianc
    yCorr      = copy(y)
    β          = zeros(nFixedEffects)  # sample of fixed effects
    α          = zeros(nMarkers)       # sample of marker effects
    δ          = zeros(nMarkers)       # inclusion indicator for marker effects
    π          = probFixed
    RinvSqrt   = sqrt(Rinv)

    #return values
    meanFxdEff = zeros(nFixedEffects)
    meanAlpha  = zeros(nMarkers)
    mdlFrq     = zeros(nMarkers)
    resVar     = zeros(chainLength)
    genVar     = zeros(chainLength)
    pi         = zeros(chainLength)
    scale      = zeros(chainLength)
