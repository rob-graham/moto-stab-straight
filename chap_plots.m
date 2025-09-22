function chap_plots(R)
figure; plot(R.v, real(R.lambda),'LineWidth',1.2);
xlabel('Speed Vx [m/s]'); ylabel('Re(\lambda)'); grid on; title('Real parts');

figure; plot(R.v, R.freq_hz,'LineWidth',1.2);
xlabel('Speed Vx [m/s]'); ylabel('Frequency [Hz]'); grid on; title('Frequencies');

figure; plot(R.v, R.zeta,'LineWidth',1.2);
xlabel('Speed Vx [m/s]'); ylabel('Damping ratio'); grid on; title('Damping Ratios');
end
