%% Made by Nout van den Bos

function [ynew,xnew] =  gaussianfit(corrmat,ypeak,xpeak)
    
    x  = xpeak-1:xpeak+1;
    y  = ypeak-1:ypeak+1;
    
    %calculate the sub-pixel precision peak, by estimating it with a
    %gaussian fit. x(2),y(2) are the middle coordinates.
    xnew = xpeak + (log(corrmat(y(2),x(1)))-log(corrmat(y(2),x(3))))/...
                    (2*log(corrmat(y(2),x(1)))-4*log(corrmat(y(2),x(2))) + ...
                    2*log(corrmat(y(2),x(3))));
    xnew = real(xnew);
    
    ynew = ypeak + (log(corrmat(y(1),x(2)))-log(corrmat(y(3),x(2))))/...
                    (2*log(corrmat(y(1),x(2)))-4*log(corrmat(y(2),x(2))) + ...
                    2*log(corrmat(y(3),x(2))));
    ynew = real(ynew);
    
end