figure;
hold on;
%xlim([-0.2,2*pi+0.2]);
%ylim([-0.2,1.2]);
set(gca,'color',[0.9,0.9,0.9]);
xlabel('time (years)','fontsize',60);
ylabel('yield','fontsize',60);
set(gca, 'FontSize', 30)

lw = 5;

exp_scale =12700;
decay_factor = 0.96;

decay_start = 10;

for i = 1:length(discounts(:,1))
  if i < decay_start
    exp_decay(i) = 1;
  else
    exp_decay(i) = decay_factor^(i-decay_start);
  end
end

plot(exp_decay'.*discounts(:,2),'color',[0,0.9,0.9],'linewidth',lw);
plot(discounts(:,2),'color',[0.3,0.3,0.7],'linewidth',lw);
plot(exp_scale*exp_decay','--','color',[0.5,0.5,0.5],'linewidth',lw);

legend('Decayed yields','Raw yields','Decay function');

xlim([0,40]);
