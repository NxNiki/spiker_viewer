function out = compute_osi(self, ori, rsp)

if self.osi_type==1
    a = sum(rsp(:).*sin(2*ori(:)/180*pi));
    b = sum(rsp(:).*cos(2*ori(:)/180*pi));
    out = sqrt(a^2 + b^2)/sum(rsp);
elseif self.osi_type==2
    out = abs(sum(rsp(:).*exp(1i*2*ori(:)))/sum(rsp(:)));
end

end