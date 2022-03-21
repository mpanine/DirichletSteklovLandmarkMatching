function Shape = insertCentralFElist(Shape, landmarks, radii)

for i = 1:length(landmarks)
    
    Shape = insertCentralFE(Shape, landmarks(i), radii(i) );
    
  
end



end

