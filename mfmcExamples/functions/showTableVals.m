function showTableVals(p, mfmc, mc)
    values = table;
    values.Computational_Budget = p;
    values.MFMC_Means = mfmc.means;
    values.MFMC_Std_Dev = mfmc.sigma;
    values.MFMC_Std_Dev_Analytical = mfmc.sigmaAnal;
    values.MFMC_MSE = mfmc.mse;
    values.MFMC_MSE_Analytical = mfmc.mseAnal;
    values.MC_Means = mc.means;
    values.MC_Std_Dev = mc.sigma;
    values.MC_Std_Dev_Analytical = mc.sigmaAnal;
    values.MC_MSE = mc.mse;
    values.MC_MSE_Analytical = mc.mseAnal;
    disp(values)
end

