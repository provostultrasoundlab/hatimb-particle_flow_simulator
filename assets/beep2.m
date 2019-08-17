function beep2
% Play a sine wave
res = 22050;
len = 0.1 * res;
hz = 720;
sound( sin( hz*(2*pi*(0:len)/res) ), res);