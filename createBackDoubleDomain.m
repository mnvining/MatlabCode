function B=createBackDoubleDomain(size_original_grid,vec_of_modes,cl)
    % Create dftmtx for extended period
    ZZ=cl*size_original_grid/2+size_original_grid;
    A=dftmtx(ZZ);
    % Pick specific modes, restrict
    B=A(:,vec_of_modes);