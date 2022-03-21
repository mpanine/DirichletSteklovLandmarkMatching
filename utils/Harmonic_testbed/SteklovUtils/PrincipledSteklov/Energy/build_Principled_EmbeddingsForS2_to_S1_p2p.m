function [S1FS2S1, S2FS2S1] = build_Principled_EmbeddingsForS2_to_S1_p2p(S1_basis, S1_products, S1_W, S2_basis, S2_products, S2_W, FmapS1S2, FmapS2S1, OrthoWeight, ProperWeight, BijectWeight, SteklovSettings)
%Builds the embedding for the unified basis principled approach.

switch SteklovSettings.NN_type

    %% Concatenated form -- is more principled
    case 'principled'
    
      S1FS2S1 = [OrthoWeight * S1_basis * FmapS1S2' ...  
                ProperWeight * S1_basis ...
                BijectWeight * S1_basis * FmapS2S1 ] ;
      
    
      S2FS2S1 = [OrthoWeight * S2_basis ...
                ProperWeight * S2_basis * FmapS1S2 ...
                BijectWeight * S2_basis ];
    
    %% Sum Form  -- Works, but really shouldn't -- FAST APPROACH
    case 'fast'
      S1FS2S1 = OrthoWeight * S1_basis * FmapS1S2' + ...  
                ProperWeight * S1_basis + ...
                BijectWeight * S1_basis * FmapS2S1 ;


      S2FS2S1 = OrthoWeight * S2_basis + ...
                ProperWeight * S2_basis * FmapS1S2 + ...
                BijectWeight * S2_basis;
                
end


end