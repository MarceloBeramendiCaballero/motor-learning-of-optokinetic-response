%% Project goals:%% Define constants for use in below cells
consts = constants();

%% Recreate Fig 1.C and 1.D
dt = 0.1;
[OKR, V, W] = simulate_okr(dt, 9 * consts.day, @training_5days, @okr, @dv, @dw);

figure(1)
subplot(2,1,1)
plot_okr(OKR, dt);

subplot(2,1,2)
plot_syn_weights(V,W,dt)

%% always training
dt = 0.1;
[OKR, V, W] = simulate_okr(dt, 9 * consts.day, @always_train, @okr, @dv, @dw);

close all

plot_okr(OKR, dt);
plot_syn_weights(V,W,dt)

%% random training (for funsies)
dt = 0.1;
[OKR, V, W] = simulate_okr(dt, 1 * consts.day, @random_training, @okr, @dv, @dw);

close all

plot_okr(OKR, dt);
plot_syn_weights(V,W,dt)

%% Recreate Fig 3
dt = 1;
[mass_OKR, ~, ~] = simulate_okr(dt, 9 * consts.day, @training_mass, @okr, @dv, @dw);

OKR_gain_vals = zeros([4,1]);

figure(3)
t = tiledlayout(8,3,'TileSpacing','Compact','Padding','Compact');

nexttile(1, [2,2])
mass_okr_plot = plot(mass_OKR, 'Color', colors().red, 'LineWidth', 1);
mass_okr_plot.XData = mass_okr_plot.XData / consts.day * dt;
xlabel("Training period (day)")
ylim([0.2, 0.8])
xlim([0,9])
ylabel("OKR gain")
legend(mass_okr_plot, "1 hour massed training")
OKR_gain_vals(1) = mass_OKR(end);

nexttile(7, [2,2])
mass_okr_plot = plot(mass_OKR,  '--', 'Color', colors().red);
mass_okr_plot.XData = mass_okr_plot.XData / consts.day * dt;
hold on
xlabel("Training period (day)")
ylim([0.2, 0.8])
xlim([0,9])
ylabel("OKR gain")

[spaced_OKR, ~, ~] = simulate_okr(dt, 9 * consts.day, @training_space_15minsx4hours, @okr, @dv, @dw);
OKR_gain_vals(2) = spaced_OKR(end);
spaced_okr_plot = plot(spaced_OKR, 'Color', colors().yellow, 'LineWidth', 1);
spaced_okr_plot.XData = spaced_okr_plot.XData / consts.day * dt;
legend(spaced_okr_plot, "15 min training x 4 hours x 1 day")

nexttile(13, [2,2])
mass_okr_plot = plot(mass_OKR, '--', 'Color', colors().red);
mass_okr_plot.XData = mass_okr_plot.XData / constants().day * dt;
xlabel("Training period (day)")
ylim([0.2, 0.8])
xlim([0,9])
ylabel("OKR gain")
hold on;

[spaced_OKR, ~, ~] = simulate_okr(dt, 9 * consts.day, @training_space_15minsx4days, @okr, @dv, @dw);
OKR_gain_vals(3) = spaced_OKR(end);
spaced_okr_plot = plot(spaced_OKR, 'Color', colors().green, 'LineWidth', 1);
spaced_okr_plot.XData = spaced_okr_plot.XData / consts.day * dt;
legend(spaced_okr_plot, "15 min training x 4 days")

nexttile(19, [2,2])
mass_okr_plot = plot(mass_OKR, '--', 'Color', colors().red);
mass_okr_plot.XData = mass_okr_plot.XData / constants().day * dt;
xlabel("Training period (day)")
xlim([0,9])
ylim([0.2, 0.8])
ylabel("OKR gain")
hold on;

[spaced_OKR, ~, ~] = simulate_okr(dt, 9 * consts.day, @training_space_7_5minsx8days, @okr, @dv, @dw);
OKR_gain_vals(4) = spaced_OKR(end);
spaced_okr_plot = plot(spaced_OKR, 'Color', colors().blue, 'LineWidth', 1);
spaced_okr_plot.XData = spaced_okr_plot.XData / consts.day * dt;
legend(spaced_okr_plot, "7.5 min training x 8 days")

nexttile(12, [5,1])
bars = bar(diag(OKR_gain_vals), 'stacked');
bars(1).FaceColor = colors().red;
bars(2).FaceColor = colors().yellow;
bars(3).FaceColor = colors().green;
bars(4).FaceColor = colors().blue;
ylim([0, 0.45])
ylabel('OKR gain')
title("Final OKR Gain")

%% Recreate Fig 4

dt = 0.1;
sim_length = 9 * consts.day;
Nt = sim_length / dt;

% 4A, deficient PF-LTP, OKR is constant since dw/dt = dv/dt = 0, so no need to run simulation
delta_OKR_PF_LTP_deficient = zeros(Nt, 1);

% 4B, deficient PF-LTD
[OKR_PF_LTD_deficient, v_ltd, w_ltd] = simulate_okr(dt, sim_length, @training_8days, @okr_def_ltd, @dv_def_ltd, @dw_def_ltd, 0, 1.1);
okr_0 = OKR_PF_LTD_deficient(1);
delta_OKR_PF_LTD_deficient = OKR_PF_LTD_deficient - okr_0; 

% 4C, selective inhibition of GABA_A receptors on PCs
[OKR_dep_gaba, v_gaba, w_gaba] = simulate_okr(dt, sim_length, @training_8days, @okr_dep_gaba, @dv_dep_gaba, @dw_dep_gaba, 0);
okr_0 = OKR_dep_gaba(1);
delta_OKR_dep_gaba = OKR_dep_gaba - okr_0;

figure(4)

subplot(3, 1, 1)
plot_okr(delta_OKR_PF_LTP_deficient, dt);
ylim([-0.1,0.6])
ylabel("\Delta OKR Gain")

subplot(3,1,2)
plot_okr(delta_OKR_PF_LTD_deficient, dt);
ylim([-0.1,0.6])
ylabel("\Delta OKR Gain")

subplot(3,1,3)
plot_okr(delta_OKR_dep_gaba, dt);
ylim([-0.1,0.6])
ylabel("\Delta OKR Gain")

figure(5)
plot_syn_weights(v_gaba, w_gaba, dt);

%% Extension: Defficient PF-LTD and depleted GABA receptors
dt = 0.1;
sim_length = 9 * consts.day;
Nt = sim_length / dt;

[OKR_PF_LTD_and_PF_LTP_deficient, v_ltd_ltp, w_ltd_ltp] = simulate_okr(dt, sim_length, @training_8days, @okr_def_ltd_def_ltp, @dv_def_ltd_def_ltp, @dw_def_ltd_def_ltp, 0, 1.1);
okr_0 = OKR_PF_LTD_and_PF_LTP_deficient(1);
delta_OKR_PF_LTD_and_PF_LTP_deficient = OKR_PF_LTD_and_PF_LTP_deficient - okr_0; 

figure(5)

subplot(2, 1, 1)
plot_okr(delta_OKR_PF_LTD_and_PF_LTP_deficient, dt);
ylim([-0.1,0.6])
ylabel("\Delta OKR Gain")

subplot(2,1,2)
plot_syn_weights(v_ltd_ltp,w_ltd_ltp,dt)

%% Training paradigms for LTD-deficiency
dt = 1;
[mass_OKR, ~, ~] = simulate_okr(dt, 2 * consts.day, @training_entire_hour_together, @okr_def_ltd, @dv_def_ltd, @dw_def_ltd, 0, 1.1);

OKR_gain_vals = zeros([4,1]);

figure()
t = tiledlayout(8,3,'TileSpacing','Compact','Padding','Compact');

nexttile(1, [2,2])
mass_okr_plot = plot(mass_OKR, 'Color', colors().red, 'LineWidth', 1);
mass_okr_plot.XData = mass_okr_plot.XData / consts.day * dt;
xlabel("Training period (day)")
ylim([-0.1, 0.4])
xlim([0,2])
ylabel("OKR gain")
OKR_gain_vals(1) = mass_OKR(end);
legend(mass_okr_plot, "1 hour massed training", "Location", "northwest")

nexttile(7, [2,2])
mass_okr_plot = plot(mass_OKR,  '--', 'Color', colors().red);
mass_okr_plot.XData = mass_okr_plot.XData / consts.day * dt;
hold on
xlabel("Training period (day)")
ylim([-0.1, 0.4])
xlim([0,2])
ylabel("OKR gain")

[spaced_OKR, ~, ~] = simulate_okr(dt, 2 * consts.day, @training_20min_every_8_hours, @okr_def_ltd, @dv_def_ltd, @dw_def_ltd, 0, 1.1);
OKR_gain_vals(2) = spaced_OKR(end);
spaced_okr_plot = plot(spaced_OKR, 'Color', colors().yellow, 'LineWidth', 1);
spaced_okr_plot.XData = spaced_okr_plot.XData / consts.day * dt;

legend(spaced_okr_plot, "20 min training x 3 sessions", "Location", "northwest")

nexttile(13, [2,2])
mass_okr_plot = plot(mass_OKR, '--', 'Color', colors().red);
mass_okr_plot.XData = mass_okr_plot.XData / constants().day * dt;
xlabel("Training period (day)")
ylim([-0.1, 0.4])
xlim([0,2])
ylabel("OKR gain")

hold on;

[spaced_OKR, ~, ~] = simulate_okr(dt, 2 * consts.day, @training_15min_every_6_hours, @okr_def_ltd, @dv_def_ltd, @dw_def_ltd, 0, 1.1);
OKR_gain_vals(3) = spaced_OKR(end);
spaced_okr_plot = plot(spaced_OKR, 'Color', colors().green, 'LineWidth', 1);
spaced_okr_plot.XData = spaced_okr_plot.XData / consts.day * dt;
legend(spaced_okr_plot, "15 min training x 4 sessions", "Location", "northwest")

nexttile(19, [2,2])
mass_okr_plot = plot(mass_OKR, '--', 'Color', colors().red);
mass_okr_plot.XData = mass_okr_plot.XData / constants().day * dt;
xlabel("Training period (day)")
ylim([-0.1, 0.4])
xlim([0,2])
ylabel("OKR gain")
hold on;

[spaced_OKR, ~, ~] = simulate_okr(dt, 2 * consts.day, @training_10min_every_4_hours, @okr_def_ltd, @dv_def_ltd, @dw_def_ltd, 0, 1.1);
OKR_gain_vals(4) = spaced_OKR(end);
spaced_okr_plot = plot(spaced_OKR, 'Color', colors().blue, 'LineWidth', 1);
spaced_okr_plot.XData = spaced_okr_plot.XData / consts.day * dt;
legend(spaced_okr_plot, "10 min training x 6 sessions", "Location", "northwest")

nexttile(12, [5,1])
bars = bar(diag(OKR_gain_vals), 'stacked');
bars(1).FaceColor = colors().red;
bars(2).FaceColor = colors().yellow;
bars(3).FaceColor = colors().green;
bars(4).FaceColor = colors().blue;
ylabel('OKR gain')
ylim([0, 0.3])


%% Effect of time constants
dt = 1;

tau_learns = 1 * consts.min : 2 * consts.min: 40 * consts.min;
tau_recovs = 1 * consts.hour : 0.25 * consts.hour : 2.5 * consts.hour;
tau_vs = 1 * consts.hour : 0.5 * consts.hour : 10 * consts.hour;

final_OKRs_for_tau_learns = zeros(1,length(tau_learns));
final_OKRs_for_tau_recovs = zeros(1,length(tau_recovs));
final_OKRs_for_tau_vs = zeros(1,length(tau_vs));

for i = 1:length(tau_learns)
    [OKR, V, W] = simulate_okr2(dt, 9 * consts.day, @training_5days, @okr, @dv_2, @dw_2, 1,1, tau_learns(i), consts.tau_recov, consts.tau_v);
    final_OKRs_for_tau_learns(i) = OKR(end);
end

for i = 1:length(tau_recovs)
    [OKR, V, W] = simulate_okr2(dt, 9 * consts.day, @training_5days, @okr, @dv_2, @dw_2, 1,1, consts.tau_learn, tau_recovs(i), consts.tau_v);
    final_OKRs_for_tau_recovs(i) = OKR(end);
end

for i = 1:length(tau_vs)
    [OKR, V, W] = simulate_okr2(dt, 9 * consts.day, @training_5days, @okr, @dv_2, @dw_2, 1,1, consts.tau_learn, consts.tau_recov, tau_vs(i));
    final_OKRs_for_tau_vs(i) = OKR(end);
end

disp("plotting")

figure

subplot(1, 3, 1)
plot(tau_learns, final_OKRs_for_tau_learns)
xlabel('\tau_{learn} (minutes)')
ylabel('Final OKR gain')
title('\tau_{learn} effect on final OKR gain')

subplot(1, 3, 2)
plot(tau_recovs/60, final_OKRs_for_tau_recovs)
xlabel('\tau_{recov} (hours)')
ylabel('Final OKR gain')
title('\tau_{recov} effect on final OKR gain')

subplot(1, 3, 3)
plot(tau_vs/60, final_OKRs_for_tau_vs)
xlabel('\tau_v (hours)')
ylabel('Final OKR gain')
title('\tau_v effect on final OKR gain')

%% Showing why tau_v has no effect on final OKR gain
dt = 0.1;

tau_vs = [1 * consts.hour,  10 * consts.hour]

[OKR1, V1, W1] = simulate_okr2(dt, 9 * consts.day, @training_5days, @okr, @dv_2, @dw_2, 1,1, consts.tau_learn, consts.tau_recov, tau_vs(1));
[OKR2, V2, W2] = simulate_okr2(dt, 9 * consts.day, @training_5days, @okr, @dv_2, @dw_2, 1,1, consts.tau_learn, consts.tau_recov, tau_vs(2));


figure

subplot(4, 1, 1)
plot_okr2(OKR1, dt)
legend('tau_v = 1 hour')

subplot(4, 1, 2)
plot_syn_weights2(V1, W1, dt)

subplot(4, 1, 3)
plot_okr2(OKR2, dt)
legend('tau_v = 1000 hours')

subplot(4, 1, 4)
plot_syn_weights2(V2, W2, dt)


%% extension trial
dt = 1;
[OKR_ex, V_ex, W_ex] = simulate_okr(dt, 9 * consts.day, @training_5days, @okr, @dv_funky, @dw);

close all

figure()
plot_okr(OKR_ex, dt)

figure()
plot_syn_weights(V_ex, W_ex, dt)


%% define training paradigms for extension




%% define parameters

function C = constants()
C.min = 1;
C.hour = 60 * C.min;
C.day = 24 * C.hour;
C.tau_learn = 20 * C.min;
C.tau_recov = 2.5 * C.hour;
C.tau_v = 5.5 * C.hour;
C.c_okr = 0.3;
C.c_okr_ltd = 1.0;
C.g_okr = 0.3;
C.g_okr_ltd = 0.3;
C.w_mli = 1.0;
C.w_0 = 1.0; 
C.w_0_ltd = 1.1;
C.v_0 = 1.0;
C.c_compensate = 1.0;
end

function colors = colors()
    colors.red = "#A2142F";
    colors.yellow = "#EDB120";
    colors.green = "#77AC30";
    colors.blue = "#0072BD";
end

%% plotting functions
function plot_okr(OKR, dt)

okr_plot = plot(OKR, 'LineWidth', 1);

okr_plot.XData = okr_plot.XData / constants().day * dt;
xlabel("Training period (day)")

ylim([0.2, 0.8])
ylabel("OKR gain")

end

function plot_syn_weights(V, W, dt)

day = constants().day;

v_plot = plot(V, 'LineWidth', 1);
hold on
w_plot = plot(W, 'LineWidth', 1);

v_plot.XData = v_plot.XData / day * dt;
w_plot.XData = w_plot.XData / day * dt;
xlabel("Training period (day)")

ylim([0, 2])
ylabel("Synaptic weight")

legend("MF-VN synaptic weight","PF-PC synaptic weight","Location","southeast")

end

%% define normal update functions

function val = okr(v, w)
    consts = constants();
    oth = v;
    oth = oth - w;
    oth = oth + consts.w_mli;
    val = consts.g_okr * (oth);
end

function val = dv(w, ~)
    consts = constants();
    val = (1.0/consts.tau_v) * (-w + consts.w_mli);
end

function val = dw(training, w)
    consts = constants();
    if training
        val = (1.0/consts.tau_learn)*(-w + consts.w_0 - consts.c_okr);
    else
        val = (1.0/consts.tau_recov)*(-w + consts.w_0);
    end
end

%% define update functions for deficient PF-LTP (Fig 4.A)
function val = okr_def_ltp(v, ~)
    val = constants().g_okr * v;
end

function val = dw_def_ltp(~, ~)
    val = 0;
end

function val = dv_def_ltp(~, ~)
    val = 0;
end

%% define update functions for deficient PF-LTD (Fig 4.B)
    function val = okr_def_ltd(v, w)
    consts = constants();
    tmp = v - w + consts.w_mli;
    if tmp >= 0
        val = tmp * consts.g_okr_ltd; % g_okr was set to one in this scenario? says so in paper but code says set c_okr to 1, not g_okr
    else
        val = 0;
    end
end

function val = dw_def_ltd(training, w)
    consts = constants();
    if training
        val = (1.0/consts.tau_learn)*(-w + consts.w_0_ltd - consts.c_okr_ltd);
    else
        val = (1.0/consts.tau_recov)*(-w + consts.w_0_ltd);
    end
end

function val = dv_def_ltd(w, v)
    consts = constants();
    tmp = -w + consts.w_mli;
    if tmp + v >= 0
        val = tmp / consts.tau_v;
    else
        val = 0;
    end
end

%% define update functions for selective depletion of GABA receptors (Fig 4.C / Wulff)
function val = okr_dep_gaba(v, w)
    consts = constants();
    oth = v;
    oth = oth - w;
    oth = oth * consts.g_okr; 
    val = oth + consts.c_compensate;
end

function val = dv_dep_gaba(w, ~)
    val = (1 / constants().tau_v) * -w;
end

function val = dw_dep_gaba(training, w)
    val = dw(training, w);
end

%% Extension: PF-LTD and PF-LTP defficiencies. Same equations as for only PF-LTD defficiency, but with wMLI = 0 (like in PF-LTP defficiency)

function val = okr_def_ltd_def_ltp(v, w)
    consts = constants();
        c_compensate = 1;
    tmp = v - w + 0;
    if tmp >= 0
        val = tmp * consts.g_okr_ltd + c_compensate; % g_okr was set to one in this scenario? says so in paper but code says set c_okr to 1, not g_okr
    else
        val = 0;
    end
end

function val = dw_def_ltd_def_ltp(training, w)
    consts = constants();
    if training
        val = (1.0/consts.tau_learn)*(-w + consts.w_0_ltd - consts.c_okr_ltd);
    else
        val = (1.0/consts.tau_recov)*(-w + consts.w_0_ltd);
    end
end

function val = dv_def_ltd_def_ltp(w, v)
    consts = constants();
    tmp = -w + 0;
    if tmp + v >= 0
        val = tmp / consts.tau_v;
    else
        val = 0;
    end
end
%% training paradigms for extension
function val = dv_funky(w, v)
    val = (1/constants().tau_v) * (-w - v / 10000);
end

function bool = training_entire_hour_together(t)

    consts = constants();
    day = consts.day;
    hour = consts.hour;
    min = consts.min;
    bool = (1*day + 22*hour + 0*min <= t && t <= 1*day + 23*hour + 0*min);
    
end


function bool = training_10min_every_4_hours(t)

    consts = constants();
    day = consts.day;
    hour = consts.hour;
    min = consts.min;
    bool = ((1*day + 2*hour +  50*min <= t && t <= 1*day + 3*hour + 0*min) ...
        || (1*day + 6*hour + 50*min <= t && t <= 1*day + 7*hour + 0*min) ...
        || (1*day + 10*hour + 50*min <= t && t <= 1*day + 11*hour + 0*min) ...
        || (1*day + 14*hour + 50*min <= t && t <= 1*day + 15*hour + 0*min)...
        || (1*day + 18*hour + 50*min <= t && t <= 1*day + 19*hour + 0*min)...
        || (1*day + 22*hour + 50*min <= t && t <= 1*day + 23*hour + 0*min));
end

function bool = training_15min_every_6_hours(t)

    consts = constants();
    day = consts.day;
    hour = consts.hour;
    min = consts.min;
    bool = ((1*day + 4*hour +  45*min <= t && t <= 1*day + 5*hour + 0*min) ...
        || (1*day + 10*hour + 45*min <= t && t <= 1*day + 11*hour + 0*min) ...
        || (1*day + 16*hour + 45*min <= t && t <= 1*day + 17*hour + 0*min) ...
        || (1*day + 22*hour + 45*min <= t && t <= 1*day + 23*hour + 0*min));
end

function bool = training_20min_every_8_hours(t)

    consts = constants();
    day = consts.day;
    hour = consts.hour;
    min = consts.min;
    bool = ((1*day + 6*hour +  40*min <= t && t <= 1*day + 7*hour + 0*min) ...
        || (1*day + 14*hour + 40*min <= t && t <= 1*day + 15*hour + 0*min) ...
        || (1*day + 22*hour + 40*min <= t && t <= 1*day + 23*hour + 0*min));
end

%% define training paradigms
% these functions define when training happens
function bool = training_5days(t)
    consts = constants();
    if (1*consts.day <= t && t < 6*consts.day && mod(t, consts.day) < 1*consts.hour)
        bool = 1;
    else
        bool = 0;
    end
end

function bool = training_8days(t)
    consts = constants();
    if (1*consts.day <= t && t < 9*consts.day && mod(t, consts.day) < 1*consts.hour)
        bool = 1;
    else
        bool = 0;
    end
end

function bool = always_train(~)
    bool = 1;
end

function bool = random_training(~)
    bool = rand() < 0.1;
end

function bool = training_mass(t)
    consts = constants();
    bool = (1*consts.day <= t && t < 2*consts.day && mod(t,consts.day) < 1*consts.hour);
end

function bool = training_space_15minsx4days(t)
    consts = constants();
    bool = (1*consts.day <= t && t < 5*consts.day && mod(t,consts.day) < 15*consts.min);
end

function bool = training_space_15minsx4hours(t)
    consts = constants();
    day = consts.day;
    hour = consts.hour;
    min = consts.min;
    bool = ((1*day + 0*hour +  0*min <= t && t <= 1*day + 0*hour + 15*min) ...
        || (1*day + 1*hour + 15*min <= t && t <= 1*day + 1*hour + 30*min) ...
        || (1*day + 2*hour + 30*min <= t && t <= 1*day + 2*hour + 45*min) ...
        || (1*day + 3*hour + 45*min <= t && t <= 1*day + 3*hour + 60*min));
end

function bool = training_space_15minsx2hours(t)
    consts = constants();
    day = consts.day;
    min = consts.min;
    bool = ((1*day +  0*min +  0*min <= t && t <= 1*day +  0*min + 15*min) ...
        || (1*day + 30*min + 15*min <= t && t <= 1*day + 30*min + 30*min) ...
        || (1*day + 60*min + 30*min <= t && t <= 1*day + 60*min + 45*min) ...
        || (1*day + 90*min + 45*min <= t && t <= 1*day + 90*min + 60*min));
end

function bool = training_space_7_5minsx8days(t)
    consts = constants();
    day = consts.day;
    min = consts.min;
    bool = (1*day <= t && t < 9*day && mod(t,day) < 7.5*min);
end

%% define simulation loop
function [OKR, V, W] = simulate_okr(dt, sim_length, is_training, okr_func, dv_func, dw_func, v0, w0)

consts = constants();

Nt = sim_length / dt;

OKR = zeros([Nt,1]);
V = zeros([Nt,1]);
W = zeros([Nt,1]);

V(1) = consts.v_0;
if nargin > 6
    V(1) = v0;
end

W(1) = consts.w_0;
if nargin > 7
    W(1) = w0;
end
OKR(1) = okr_func(V(1), W(1));

for j = 2:Nt

    V(j) = max(0, V(j-1) + dt * dv_func(W(j-1), V(j-1)));
    W(j) = max(0, W(j-1) + dt * dw_func(is_training((j-1) * dt), W(j-1)));
    OKR(j) = okr_func(V(j), W(j));

end

end


%% Defining functions for extension evaluating the effect of time constants

function val = dv_2(w, ~, tau_v)
    consts = constants();
    val = (1.0/tau_v) * (-w + consts.w_mli);
end

function val = dw_2(training, w, tau_learn, tau_recov)
    consts = constants();
    if training
        val = (1.0/tau_learn)*(-w + consts.w_0 - consts.c_okr);
    else
        val = (1.0/tau_recov)*(-w + consts.w_0);
    end
end

%% define simulate_okr for extension. It takes the set of constants as an input
function [OKR, V, W] = simulate_okr2(dt, sim_length, is_training, okr_func, dv_func, dw_func, v0, w0, tau_learn, tau_recov, tau_v)

consts = constants();

Nt = sim_length / dt;

OKR = zeros([Nt,1]);
V = zeros([Nt,1]);
W = zeros([Nt,1]);

V(1) = consts.v_0;
if nargin > 6
    V(1) = v0;
end

W(1) = consts.w_0;
if nargin > 7
    W(1) = w0;
end
OKR(1) = okr_func(V(1), W(1));

for j = 2:Nt

    V(j) = max(0, V(j-1) + dt * dv_func(W(j-1), V(j-1), tau_v));
    W(j) = max(0, W(j-1) + dt * dw_func(is_training((j-1) * dt), W(j-1), tau_learn, tau_recov));
    OKR(j) = okr_func(V(j), W(j));

end

end

%% Plotting functions for changing tau_v plots
function plot_okr2(OKR, dt)

okr_plot = plot(OKR, 'LineWidth', 1);

okr_plot.XData = okr_plot.XData / constants().day * dt;
xlabel("Training period (day)")

ylim([0.2, 0.95])
ylabel("OKR gain")

end

function plot_syn_weights2(V, W, dt)

day = constants().day;

v_plot = plot(V, 'LineWidth', 1);
hold on
w_plot = plot(W, 'LineWidth', 1);

v_plot.XData = v_plot.XData / day * dt;
w_plot.XData = w_plot.XData / day * dt;
xlabel("Training period (day)")

ylim([0, 3])
ylabel("Synaptic weight")

legend("MF-VN synaptic weight","PF-PC synaptic weight","Location","southeast")

end