function [W]=dBm_W(dBm)
    W=10^((dBm-30)/10);
end
