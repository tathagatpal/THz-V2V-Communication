function [g1, g2] = chan_gen(fc,N,k)
c = 3*10^8;
lambda = c/(fc);

for i = 1:N
    g1_los = sqrt(0.5)*(randn+1i*randn)*(1/sqrt(2));
    g1_nlos_rician = sqrt(0.5)*(randn+1i*randn)*(1/sqrt(2));
    g1(i) = sqrt(k/(k+1))*g1_nlos_rician + ...
        sqrt(k/(k+1))*g1_los;
end

for i = 1:N
    g2_los = sqrt(0.5)*(randn+1i*randn)*(1/sqrt(2));
    g2_nlos_rician = sqrt(0.5)*(randn+1i*randn)*(1/sqrt(2));
    g2(i) = sqrt(k/(k+1))*g2_nlos_rician  + ...
        sqrt(k/(k+1))*g2_los;
end

end