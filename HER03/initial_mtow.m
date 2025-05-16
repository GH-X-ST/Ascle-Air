c = 0.25; % specific fuel consumption kg/kWh
P = 500; %power kW
dto = 0.05; % 3 min take off
dc = 0.417; % 25 min cruise
dl = 0.05; % 3 min land
we_by_w0 = 0.5; %from raymer empty to total weight
wpayload = 1200; % estimate kg
numtrips = 4;

w0 = 3500; %kg initial guess

% take off -> cruise -> land -> take off -> cruise -> land

for i =1:50

    wto1 = 1 - (c*dto)/(w0/P);
    wc1 = 1 - (c*dc)/(wto1*w0/P);
    wl1 = 1 - (c*dl)/(wc1*w0/P);
    wto2 = 1 - (c*dto)/(wl1*w0/P);
    wc2 = 1 - (c*dc)/(wto2*w0/P);
    wl2 = 1 - (c*dl)/(wc2*w0/P);
    
    nosafetywfbyw0 = (wto1*wc1*wl1*wto2*wc2*wl2) ^ numtrips;
    wsafety = 1 - (c*0.3333)/(wl2*w0/P);
    wf_by_w0 = 1 - (wto1*wc1*wl1*wto2*wc2*wl2) ^ numtrips *wsafety;
  


    w0 = (wpayload)/(1-wf_by_w0-we_by_w0);
end
