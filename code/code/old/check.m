figure(1234)
for channel=1:16;
subplot(4,4,channel);
imagesc(reshape(Phi(channel, :, :), M, w))
az = 0;
el = 90;
view(az, el);
colormap('gray');
end