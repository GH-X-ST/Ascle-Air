framp = 1;
nz = 3.5*1.5; % ultimate load factor (*1.5?)
Sbody = 517; %ft squared this is a bad estimate maybe check cad - its total wetted area so of the whole body
l = 13.3; % length just a guess
comb_no_fuselage = 3119;
fuselage_guess = 450;

for i =1:10
    wbasic = 5.896*framp*((comb_no_fuselage+fuselage_guess)/1000)^0.4908*nz^0.1323*Sbody^0.2544*l^0.61;
    fuselage_guess = wbasic
end
