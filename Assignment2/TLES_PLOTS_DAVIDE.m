
%% Plot from TLEs
% figure(4)
%
% [DATES, KEP]=readTLE('debris5.tle', 0);
%
% KEP(:,3:5) = unwrap(KEP(:,3:5),[],1);
% KEP(5:end,6) = unwrap(KEP(5:end,6),[],1);
%
% KEP(:,3:6) = rad2deg(KEP(:,3:6));
% % Semi-major axis
% subplot(2,3,1)
%
% plot(datenum(DATES),KEP(:,1));
% legend('TLEs')
% grid on
% xlabel('${time [T]}$','Interpreter', 'latex','Fontsize', 14)
% ylabel('$\mathbf{a [Km]}$','Interpreter', 'latex','Fontsize', 14)
%
% % Eccentricity
% subplot(2,3,2)
% plot(datenum(DATES),KEP(:,2));
% legend('TLEs')
% grid on
% xlabel('${time [T]}$','Interpreter', 'latex','Fontsize', 14)
% ylabel('$\mathbf{e [-]}$','Interpreter', 'latex','Fontsize', 14)
%
%
% % inclination
% subplot(2,3,3)
% plot(datenum(DATES),KEP(:,3));
% legend('TLEs')
% grid on
% xlabel('${time [T]}$','Interpreter', 'latex','Fontsize', 14)
% ylabel('$\mathbf{i [deg]}$','Interpreter', 'latex','Fontsize', 14)
%
% % RAAN
% subplot(2,3,4)
% plot(datenum(DATES),KEP(:,4));
% legend('TLEs')
% grid on
% xlabel('${time [T]}$','Interpreter', 'latex','Fontsize', 14)
% ylabel('$\mathbf{\Omega  [deg]}$','Interpreter', 'latex','Fontsize', 14)
%
% % omega
% subplot(2,3,5)
% plot(datenum(DATES),KEP(:,5));
% legend('TLEs')
% grid on
% xlabel('${time [T]}$','Interpreter', 'latex','Fontsize', 14)
% ylabel('$\mathbf{\omega  [deg]}$','Interpreter', 'latex','Fontsize', 14)
%
%
% % f
% subplot(2,3,6)
% plot(datenum(DATES),KEP(:,6));
% legend('TLEs')
% grid on
% xlabel('${time [T]}$','Interpreter', 'latex','Fontsize', 14)
% ylabel('$\mathbf{f  [deg]}$','Interpreter', 'latex','Fontsize', 14)
