function generateMassMatrix(FEMmesh)
    M = Mass2DMatrix(FEMmesh)
    #Mlu = lu(M)
    return M
end

function generateMassMatrixWScalar(FEMmesh,f)
    M = Mass2DMatrix(FEMmesh,f)
    #Mlu = lu(M)
    return M
end
