function y=my_smooth(data, z, Wsize)

window=ones(Wsize,1)./Wsize;
y=filtfilt(window,1,data);