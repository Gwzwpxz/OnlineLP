function [data] = gendata(m, n, type, tightness, sps, bord)
% Data generator for online algorithm test
% [data] = gendata(m, n, type, tightness, sps, bord)
% Only two types of data will be considered during trace of the test:
% 1. Chu & Beasley's MKP Benchmark
% 2. LSY's data generation method

if type == 1
    
    if sps < 1
        data.A = round(sprand(m, n, sps), 2) * 10;
    else
        data.A = randi(1000, m, n) / 100;
    end % End if
    data.b = full(tightness * sum(data.A, 2) / (n^(1 - bord)));
    data.c = full(sum(data.A, 1)'/m + rand(n, 1) * 5);

else
    
    data.A = randn(m, n) + 1;
    data.b = n^bord * (rand(m, 1) + 1) / 3;
    data.c = sum(data.A, 1)'/m - rand(n, 1) * m;
    
end % End if

end % End function
