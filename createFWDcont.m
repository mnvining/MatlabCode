function C=createFWDcont(size_original_grid,vec_of_modes,cl,ZZ)
    %ZZ=cl*size_original_grid/2+size_original_grid;

    % Create dftmtx for extended period
    A=dftmtx(ZZ);
    A=A';
    A=A/sqrt(ZZ);
    % Pick specific modes, restrict
    C=A(1:size_original_grid,vec_of_modes);