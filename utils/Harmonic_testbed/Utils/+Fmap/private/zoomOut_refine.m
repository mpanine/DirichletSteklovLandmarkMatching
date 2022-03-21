function refined = zoomOut_refine( BASIS_Src, BASIS_Tar, pF, startSize, step_size )
C21in = fmap_encode_nn(BASIS_Src.basis, BASIS_Tar.basis, pF, startSize, startSize);

dim_out = size(BASIS_Src.basis, 2);

Cmaps = fmap_upsample_keep_all(BASIS_Src.basis, BASIS_Tar.basis, C21in, dim_out, step_size);
Cout = Cmaps{end};

refined = fmap_decode_nn(BASIS_Src.basis, BASIS_Tar.basis, Cout);
end

function all_Cs = fmap_upsample_keep_all(X1, X2, Cin, dim_out, step)
    C = Cin;

    all_Cs = {};
    all_Cs{1} = C;
    for dim=size(C,2)+1:step:dim_out
        C = fmap_encode_nn(X1,X2, fmap_decode_nn(X1,X2,C), dim, dim);
        all_Cs{end+1} = C;
    end
end


function C12 = fmap_encode_nn(X1, X2, nn, dim_out1, dim_out2)
        X1out = X1(:,1:dim_out1);
        X2out = X2(:,1:dim_out2);
        
        C12 = X2out\X1out(nn,:);
end

function nn = fmap_decode_nn(X1, X2, C12)
    dim1 = size(C12,1);
    dim2 = size(C12,2);
    nn = knnsearch(X1(:,1:dim2)*C12', X2(:,1:dim1));
end