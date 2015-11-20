#BayesC
function sample_marker_effects!(method::BayesC,marker_matrix,estimated_parameters,data::Data)

    xArray  = marker_matrix.xArray
    XpRinvX = marker_matrix.XpRinvX
    α       = estimated_parameters.α

    logPi         = log(estimated_parameters.π)
    logPiComp     = log(1.0-estimated_parameters.π)
    logVarEffects = log(varEffects)
    logDelta0     = logPi
    invVarRes     = 1.0/estimated_parameters.vare
    invVarEffects = 1.0/estimated_parameters.varEffects
    nLoci = 0

    for j=1:nMarkers
        x = xArray[j]
        rhs = (dot(x.*Rinv,yCorr) + XpRinvX[j]*α[j])*invVarRes
        lhs = XpRinvX[j]*invVarRes + invVarEffects
        invLhs = 1.0/lhs
        gHat   = rhs*invLhs
        logDelta1  = -0.5*(log(lhs) + logVarEffects - gHat*rhs) + logPiComp
        probDelta1 = 1.0/(1.0 + exp(logDelta0 - logDelta1))
        oldAlpha = α[j]

        if(rand()<probDelta1)
            δ[j] = 1
            α[j] = gHat + randn()*sqrt(invLhs)
            BLAS.axpy!(oldAlpha-α[j],x,yCorr)
            nLoci = nLoci + 1
        else
            if (oldAlpha!=0)
                BLAS.axpy!(oldAlpha,x,yCorr)
            end
            δ[j] = 0
            α[j] = 0
        end
    end
    estimated_parameters.nLoci=nLoci
end

#BayesB
function sample_marker_effects!(method::BayesB,marker_matrix,estimated_parameters,data::Data)
    xArray  = marker_matrix.xArray
    XpRinvX = marker_matrix.XpRinvX
    α       = estimated_parameters.α
    u       = estimated_parameters.u

    logPi         = log(π)
    logPiComp     = log(1.0-π)
    logDelta0     = logPi
    invVarRes     = 1.0/vare
    invlocusEffectVar = 1.0./locusEffectVar
    loglocusEffectVar = log(locusEffectVar)
    nLoci = 0

    for j=1:nMarkers
        x = xArray[j]
        rhs = (dot(x.*Rinv,yCorr) + XpRinvX[j]*u[j])*invVarRes
        lhs = XpRinvX[j]*invVarRes + invlocusEffectVar[j]
        invLhs = 1.0/lhs
        gHat   = rhs*invLhs
        logDelta1  = -0.5*(log(lhs) + loglocusEffectVar[j] - gHat*rhs) + logPiComp
        probDelta1 = 1.0/(1.0 + exp(logDelta0 - logDelta1))
        oldu = u[j]

        if(rand()<probDelta1)
            δ[j] = 1
            α[j] = gHat + randn()*sqrt(invLhs)
            u[j] = α[j]
            BLAS.axpy!(oldu-u[j],x,yCorr)
            nLoci = nLoci + 1
        else
            if (oldu!=0)
                BLAS.axpy!(oldu,x,yCorr)
            end
            δ[j] = 0
            α[j] = randn()*sqrt(locusEffectVar[j])
            u[j] = 0
        end
    end
    estimated_parameters.nLoci=nLoci
end