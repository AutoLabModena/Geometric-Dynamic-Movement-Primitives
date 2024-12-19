%
% Written by Giovanni Braglia and Davide Tebaldi, 2023
% University of Modena and Reggio Emilia
% website: https://www.automatica.unimore.it/
%
%---------------------------------------------------------------------
function [tn,sn,xn]= SpatialSampling(t,x,delta)
%---------------------------------------------------------------------
% Returns equally spaced ('delta') points on trajectory 'x'

    [nrow,ncol]=size(x);
    if ncol>nrow
        x=x';
    end
    %
    for jj=1:size(x,2)
        eval(['xn' num2str(jj) '=[x(1,' num2str(jj) ')];']);
    end
    %
    sn=[t(1)];
    tn=[t(1)];
    %
    xcurr=x(1,:);
    for ii=1:size(x,1)-1
        %
        if norm(x(ii+1,:)-xcurr)>delta
            nr_steps=floor(norm(x(ii+1,:)-xcurr)/delta);
            for pp=1:nr_steps
                sn=[sn; sn(end)+delta];
                for jj=1:size(x,2)
                    xincr=xcurr(jj)+delta*(x(ii+1,jj)-xcurr(jj))/norm(x(ii+1,:)-xcurr);
                    eval(['xn' num2str(jj) '=[xn' num2str(jj) '; ' num2str(xincr) '];']);
                end
                xcurr=x(1,:)*0;
                for jj=1:size(x,2)
                    eval(['xcurr(' num2str(jj) ')=xn' num2str(jj) '(end);']);
                end
                nDX=norm(x(ii+1,:)-x(ii,:));
                DT=t(ii+1)-t(ii);
                K=nDX/DT;
                tn=[tn; (norm(xcurr-x(ii,:))+K*t(ii))/K];
            end
        end
    end
    xn=zeros(size(xn1,1),size(x,2));
    for jj=1:size(x,2)
        eval(['xn(:,' num2str(jj) ')=xn' num2str(jj) ';']);
    end

end
