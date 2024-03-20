function PerpMag = TransverseMagnetizationMagnitude_6vec(M_t)
% calculate perpendicular magnetization from a group of isochromats

avgx = mean(M_t(1,:));
avgy = mean(M_t(3,:));
PerpMag = sqrt(avgx^2 +avgy^2);

% Added for Magnetization prepared sequences, to get sign of magnetization
theta = atan2(avgy, avgx);
if theta < 0
    PerpMag = PerpMag * -1;
end
    
