function Q = extract_Q(filename)

load(filename);

lambda = scandata.wavelength;
Tlog = scandata.power(:,2)';

%Fit 1550nm peak to lorentzian

% Find index of 1550nm
[~,ind1550] = min(abs(lambda-1545e-9));

% Find peak closest to 1550nm
[~,peaks_index] = findpeaks(Tlog,'MinPeakProminence',20);
[~, indcent] = min(abs(lambda(peaks_index)-lambda(ind1550)));

indcent = peaks_index(indcent);

%Fit to data lorentzian
num = 40; 
tempT = Tlog(indcent-num:indcent+num);
templ = lambda(indcent-num:indcent+num);

%Set wavelength input to lorentz fit function to 1e12x
x = templ*1e12;
y = 10.^(tempT/10);
[~,P,~,~] = lorentzfit(x,y,[],[],'3');

figure; hold on;
scatter(x/1e6,10*log10(y));
templc = 1e12*linspace(lambda(indcent-num),lambda(indcent+num),10000);
lorentzT = 10*log10(P(1)./((templc - P(2)).^2 + P(3)));
plot(templc/1e6, 10*log10(P(1)./((templc - P(2)).^2 + P(3))));

%Analyze Lorentzian Fit data
lambda_fit = templc/1e12;
[max_trans, indQ] = max(lorentzT);

%Find 3dB frequency (BW)
[~, ind3dB] = min(abs(lorentzT-(max_trans-3)));

BW = 2*abs(lambda_fit(indQ)-lambda_fit(ind3dB));
Q = lambda_fit(indQ)/BW;
end