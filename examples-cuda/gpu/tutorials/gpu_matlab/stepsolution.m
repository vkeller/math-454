function [vv,vvold] = stepsolution(vv,vvold,ii,N,dt,W1T,W2T,W3T,W4T,...
                                        WuxxT1,WuxxT2,WuyyT1,WuyyT2)
    V = [vv(ii,:) vv(ii,N:-1:2)];
    U = real(fft(V.')).';
    
    W1test = (U.*W1T).';
    W2test = (U.*W2T).';
    W1 = (real(ifft(W1test))).';
    W2 = (real(ifft(W2test))).';
    
    % Calculating 2nd derivative in x
    uxx(ii,ii) = W2(:,ii).* WuxxT1 - W1(:,ii).*WuxxT2;
    uxx([1,N+1],[1,N+1]) = 0;
    
    V = [vv(:,ii); vv((N:-1:2),ii)];
    U = real(fft(V));
    
    W1 = real(ifft(U.*W3T));
    W2 = real(ifft(U.*W4T));
    
    % Calculating 2nd derivative in y
    uyy(ii,ii) = W2(ii,:).* WuyyT1 - W1(ii,:).*WuyyT2;
    uyy([1,N+1],[1,N+1]) = 0;
    
    % Computing new value using 2nd order central finite difference in
    % time
    vvnew = 2*vv - vvold + dt*dt*(uxx+uyy);
    vvold = vv; vv = vvnew;

end