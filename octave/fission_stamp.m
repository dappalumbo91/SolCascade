% Fission stamp. Run in Octave 11.3.0.
% Recompute each installed closing from the pinned seeds.
% A unique close meets the helium quotient match, and the next dressing stays wider.
% The four held pairs keep both routes.
% The four-hour year share is an installed close. Four factors stay wider.
% The file ledger quantum is an installed close.
% The comparison column's printed quotients are installed closes and stay unrun.
% The Red Book figures off that column and the CANDU figures are installed closes.
% The catalog ceiling, the 6.5 percent page, and the 60 GWd/t stack stay unrun.
% Domain scalars stay out of the joules.
% Not part of lean_lib FSOT.

function fission_stamp(outname)
if nargin < 1 || isempty(outname)
  outname = 'data/fission_stamp.txt';
end
S.e = 2.7182818284590451;
S.pi = 3.1415926535897931;
S.phi = 1.6180339887498949;
S.gamma = 0.57721566490153287;
S.G = 0.91596559417721901;
S.alpha = 0.0072973525649202339;
helium = helium_match(S);
if abs(helium - 4.2343138989270577e-06) / helium > 1e-7
  error('the helium quotient match left the pin');
end
[base_names, base_values] = base_pieces(S);
if length(base_names) ~= 60
  error('the base piece count changed');
end
fid = fopen(outname, 'w');
if fid < 0
  error('the stamp log did not open');
end
fails = 0;
fprintf(fid, 'helium %.16e\n', helium);
fprintf(fid, 'base_pieces %d\n', length(base_names));
fails = check_close(fid, fails, S, base_names, base_values, helium, 'decay_one_second_fraction', 0.058220348783608857, 'unit', 116, false, {'(phi^-6/sqrt(G))*(1-alpha*e^-4)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'decay_two_day_limit', 0.005953680887458567, 'unit', 116, false, {'(1-(gamma^6+sqrt(G)))*(1-alpha*pi^-5)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'decay_two_day_drop', 0.79907641804468166, 'unit', 110, false, {'(1-(gamma^3/sqrt(G)))*(1+alpha*pi^-5)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'decay_one_second_limit', 0.058889212030811468, 'unit', 116, false, {'((pi^-5+phi^-6)*(1-alpha*phi^-3))*(1-alpha*pi^-4)' '(((pi^-5+phi^-6)*(1-alpha*phi^-3))*(1-alpha*pi^-4))*(1-alpha*pi^-6)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'decay_one_day_fraction', 0.0066158926487938407, 'unit', 116, false, {'((pi^-4*G^5)*(1-alpha*gamma^5))*(1-alpha*pi^-6)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'decay_one_day_limit', 0.0072838608269594068, 'unit', 116, false, {'((e^-3*phi^-4)*(1+alpha*phi^-2))*(1-alpha*pi^-5)' '(((e^-3*phi^-4)*(1+alpha*phi^-2))*(1-alpha*pi^-5))*(1-alpha*pi^-6)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'decay_two_day_fraction', 0.0052866037999663232, 'unit', 116, false, {'((pi^-5*phi)*(1-alpha*e^-4))*(1-alpha*pi^-6)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'decay_one_day_drop', 0.11363540045738192, 'unit', 116, false, {'((gamma^3*G^6)*(1+alpha*gamma^5))*(1+alpha*e^-5)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'heat_ratio', 3.0303030303030298, 'positive', 116, false, {'(pi-gamma^4)*(1-alpha*pi^-4)' '((pi-gamma^4)*(1-alpha*pi^-4))*(1-alpha*e^-6)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'removal_cap', 0.050000000000000003, 'unit', 116, false, {'(pi^-3*G^-5)*(1-alpha*phi^-6)' '((pi^-3*G^-5)*(1-alpha*phi^-6))*(1-alpha*e^-6)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'removal_duty', 0.003738440400341451, 'unit', 116, false, {'(pi^-2*gamma^6)*(1-alpha*gamma^2)' '((pi^-2*gamma^6)*(1-alpha*gamma^2))*(1+alpha*e^-5)' '(((pi^-2*gamma^6)*(1-alpha*gamma^2))*(1+alpha*e^-5))*(1-alpha*pi^-6)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'removal_one_second', 15.573519006953626, 'positive', 116, false, {'(gamma^-5-gamma^6)*(1+alpha*pi^-3)' '((gamma^-5-gamma^6)*(1+alpha*pi^-3))*(1+alpha*pi^-6)' '(((gamma^-5-gamma^6)*(1+alpha*pi^-3))*(1+alpha*pi^-6))*(1+alpha*pi^-6)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'removal_two_day', 1.4141229138710416, 'positive', 116, false, {'(G^3+G^5)*(1+alpha*phi^-5)' '((G^3+G^5)*(1+alpha*phi^-5))*(1-alpha*e^-5)' '(((G^3+G^5)*(1+alpha*phi^-5))*(1-alpha*e^-5))*(1+alpha*pi^-6)' '((((G^3+G^5)*(1+alpha*phi^-5))*(1-alpha*e^-5))*(1+alpha*pi^-6))*(1+alpha*pi^-6)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'service_service', 0.011048726486417693, 'unit', 116, false, {'(sqrt(e)*e^-5)*(1-alpha*sqrt(gamma))' '((sqrt(e)*e^-5)*(1-alpha*sqrt(gamma)))*(1+alpha*e^-4)' '(((sqrt(e)*e^-5)*(1-alpha*sqrt(gamma)))*(1+alpha*e^-4))*(1-alpha*e^-6)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'service_cooldown', 0.0069976441470346387, 'unit', 116, false, {'(pi^-4-pi^-5)*(1-alpha*pi^-4)' '((pi^-4-pi^-5)*(1-alpha*pi^-4))*(1-alpha*pi^-6)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'service_cross', 0.014246340407594619, 'unit', 116, false, {'(phi^-2*gamma^6)*(1+alpha*G^-2)' '((phi^-2*gamma^6)*(1+alpha*G^-2))*(1-alpha*gamma^6)' '(((phi^-2*gamma^6)*(1+alpha*G^-2))*(1-alpha*gamma^6))*(1+alpha*e^-6)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'discharge', 0.047356994679361296, 'unit', 116, false, {'(e^-3-e^-6)*(1+alpha*phi^-4)' '((e^-3-e^-6)*(1+alpha*phi^-4))*(1-alpha*pi^-5)' '(((e^-3-e^-6)*(1+alpha*phi^-4))*(1-alpha*pi^-5))*(1-alpha*pi^-6)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'reload_companion', 0.98864383561643832, 'positive', 116, false, {'(sqrt(gamma)*G^-3)*(1+alpha*e^-6)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'reload_separative', 7.4897119341563787, 'positive', 116, false, {'(pi^-2+e^2)*(1-alpha*pi^-4)' '((pi^-2+e^2)*(1-alpha*pi^-4))*(1-alpha*e^-6)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'reload_highly_enriched', 26.703296703296704, 'positive', 116, false, {'(phi^5+gamma^-5)*(1+alpha*pi^-3)' '((phi^5+gamma^-5)*(1+alpha*pi^-3))*(1+alpha*pi^-6)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'reload_energy_return', 951.70549450549436, 'positive', 116, false, {'(pi^6-gamma^-4)*(1-alpha*pi^-2)' '((pi^6-gamma^-4)*(1-alpha*pi^-2))*(1+alpha*pi^-5)' '(((pi^6-gamma^-4)*(1-alpha*pi^-2))*(1+alpha*pi^-5))*(1+alpha*pi^-6)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'steam_built', 0.35999999999999999, 'positive', 116, false, {'(pi^-3*phi^5)*(1+alpha*G)' '((pi^-3*phi^5)*(1+alpha*G))*(1-alpha*e^-4)' '(((pi^-3*phi^5)*(1+alpha*G))*(1-alpha*e^-4))*(1-alpha*e^-5)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'steam_carnot', 0.4464285714285714, 'positive', 116, false, {'(phi^-2+gamma^5)*(1+alpha*gamma^4)' '((phi^-2+gamma^5)*(1+alpha*gamma^4))*(1+alpha*e^-5)' '(((phi^-2+gamma^5)*(1+alpha*gamma^4))*(1+alpha*e^-5))*(1+alpha*pi^-6)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'steam_carnot_step', 1.3528138528138527, 'positive', 116, false, {'(e^-3+G^-3)*(1+alpha*gamma^3)' '((e^-3+G^-3)*(1+alpha*gamma^3))*(1-alpha*pi^-4)' '(((e^-3+G^-3)*(1+alpha*gamma^3))*(1-alpha*pi^-4))*(1-alpha*e^-6)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'margin_balance', 0.14999999999999999, 'positive', 116, false, {'(e^-3*gamma^-2)*(1+alpha*gamma)' '((e^-3*gamma^-2)*(1+alpha*gamma))*(1-alpha*phi^-6)' '(((e^-3*gamma^-2)*(1+alpha*gamma))*(1-alpha*phi^-6))*(1+alpha*pi^-6)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'margin_fragment', 0.84999999999999998, 'positive', 116, false, {'(phi^-5+sqrt(gamma))*(1+alpha*pi^-4)' '((phi^-5+sqrt(gamma))*(1+alpha*pi^-4))*(1+alpha*pi^-5)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'margin_uranium', 0.0066, 'unit', 116, false, {'(pi^-4*G^5)*(1-alpha*phi^-2)' '((pi^-4*G^5)*(1-alpha*phi^-2))*(1-alpha*pi^-4)' '(((pi^-4*G^5)*(1-alpha*phi^-2))*(1-alpha*pi^-4))*(1-alpha*e^-6)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'margin_plutonium', 0.0023, 'unit', 116, false, {'(pi^-5*G^4)*(1-alpha*pi^-4)' '((pi^-5*G^4)*(1-alpha*pi^-4))*(1-alpha*e^-6)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'margin_uranium_233', 0.27000000000000002, 'positive', 116, false, {'(phi^-5*gamma^-2)*(1-alpha*gamma^2)' '((phi^-5*gamma^-2)*(1-alpha*gamma^2))*(1+alpha*pi^-4)' '(((phi^-5*gamma^-2)*(1-alpha*gamma^2))*(1+alpha*pi^-4))*(1+alpha*pi^-6)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'bounds_ceiling', 0.45000000000000001, 'positive', 116, true, {'(sqrt(gamma)*G^6)*(1+alpha*phi^-2)' '((sqrt(gamma)*G^6)*(1+alpha*phi^-2))*(1+alpha*e^-4)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'bounds_headroom', 1.3636363636363635, 'positive', 116, false, {'(sqrt(phi)+phi^-5)*(1+alpha*phi^-4)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'bounds_fast_efficiency', 0.40000000000000002, 'positive', 116, false, {'(phi^-3*G^-6)*(1+alpha*phi^-5)' '((phi^-3*G^-6)*(1+alpha*phi^-5))*(1+alpha*pi^-5)' '(((phi^-3*G^-6)*(1+alpha*phi^-5))*(1+alpha*pi^-5))*(1+alpha*pi^-6)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'bounds_fast_ceiling', 0.5, 'positive', 116, false, {'(sqrt(pi)-sqrt(phi))*(1-alpha*gamma^4)' '((sqrt(pi)-sqrt(phi))*(1-alpha*gamma^4))*(1-alpha*e^-5)' '(((sqrt(pi)-sqrt(phi))*(1-alpha*gamma^4))*(1-alpha*e^-5))*(1-alpha*pi^-6)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'bounds_fast_headroom', 1.25, 'positive', 116, false, {'(gamma^2+G)*(1+alpha*phi^-5)' '((gamma^2+G)*(1+alpha*phi^-5))*(1+alpha*pi^-5)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'bounds_small_class', 0.29999999999999999, 'positive', 116, false, {'(phi^-3+gamma^5)*(1-alpha*gamma^5)' '((phi^-3+gamma^5)*(1-alpha*gamma^5))*(1-alpha*pi^-6)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'bounds_late_band', 0.10000000000000001, 'unit', 116, false, {'(e-phi^2)*(1-alpha*gamma^2)' '((e-phi^2)*(1-alpha*gamma^2))*(1-alpha*e^-5)' '(((e-phi^2)*(1-alpha*gamma^2))*(1-alpha*e^-5))*(1+alpha*pi^-6)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'bounds_thermal_return', 752.58154030480125, 'positive', 116, false, {'(pi^6/sqrt(phi))*(1-alpha*gamma)' '((pi^6/sqrt(phi))*(1-alpha*gamma))*(1-alpha*e^-5)' '(((pi^6/sqrt(phi))*(1-alpha*gamma))*(1-alpha*e^-5))*(1+alpha*pi^-6)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'bounds_fast_return', 896.92556081066539, 'positive', 116, false, {'(pi^6-e^4)*(1-alpha*G^-5)' '((pi^6-e^4)*(1-alpha*G^-5))*(1+alpha*gamma^5)' '(((pi^6-e^4)*(1-alpha*G^-5))*(1+alpha*gamma^5))*(1-alpha*pi^-5)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'bounds_efficiency_step', 1.2121212121212122, 'positive', 116, false, {'(e^-4+G^-2)*(1+alpha*phi^-3)' '((e^-4+G^-2)*(1+alpha*phi^-3))*(1-alpha*e^-4)' '(((e^-4+G^-2)*(1+alpha*phi^-3))*(1-alpha*e^-4))*(1-alpha*e^-6)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'stack_atomic', 0.065000000000000002, 'unit', 116, true, {'(pi^-6+gamma^5)*(1-alpha*phi^-3)' '((pi^-6+gamma^5)*(1-alpha*phi^-3))*(1-alpha*e^-5)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'stack_burnup', 1.3333333333333333, 'positive', 116, true, {'(pi^-3+G^-3)*(1-alpha*e^-4)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'stack_stack', 1.8037518037518034, 'positive', 116, true, {'(gamma^4+G^-6)*(1-alpha*gamma^6)' '((gamma^4+G^-6)*(1-alpha*gamma^6))*(1-alpha*e^-6)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'stack_higher', 0.063142659572481732, 'unit', 116, true, {'(e^-3*sqrt(phi))*(1-alpha*phi^-2)' '((e^-3*sqrt(phi))*(1-alpha*phi^-2))*(1-alpha*e^-4)' '(((e^-3*sqrt(phi))*(1-alpha*phi^-2))*(1-alpha*e^-4))*(1-alpha*e^-5)' '((((e^-3*sqrt(phi))*(1-alpha*phi^-2))*(1-alpha*e^-4))*(1-alpha*e^-5))*(1+alpha*pi^-6)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'stack_electric', 0.085420264573162513, 'unit', 116, true, {'(gamma^4*G^3)*(1+alpha*gamma^3)' '((gamma^4*G^3)*(1+alpha*gamma^3))*(1-alpha*pi^-4)' '(((gamma^4*G^3)*(1+alpha*gamma^3))*(1-alpha*pi^-4))*(1-alpha*pi^-6)' '((((gamma^4*G^3)*(1+alpha*gamma^3))*(1-alpha*pi^-4))*(1-alpha*pi^-6))*(1-alpha*pi^-6)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'stack_enrichment', 4.5, 'positive', 116, false, {'(gamma^-3-G^4)*(1+alpha*e^-2)' '((gamma^-3-G^4)*(1+alpha*e^-2))*(1-alpha*pi^-4)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'fit_weight', 0.87, 'positive', 116, false, {'(pi^-2+G^3)*(1+alpha*pi^-3)' '((pi^-2+G^3)*(1+alpha*pi^-3))*(1-alpha*e^-6)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'fit_exponent', 0.20000000000000001, 'positive', 116, false, {'(e^-2+gamma^5)*(1+alpha*phi^-2)' '((e^-2+gamma^5)*(1+alpha*phi^-2))*(1+alpha*e^-4)' '(((e^-2+gamma^5)*(1+alpha*phi^-2))*(1+alpha*e^-4))*(1+alpha*pi^-5)' '((((e^-2+gamma^5)*(1+alpha*phi^-2))*(1+alpha*e^-4))*(1+alpha*pi^-5))*(1+alpha*pi^-6)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'condenser_one_second', 0.086896384723691905, 'unit', 116, false, {'(e^-3+gamma^6)*(1+alpha*gamma^3)' '((e^-3+gamma^6)*(1+alpha*gamma^3))*(1+alpha*pi^-5)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'condenser_one_day', 0.0098744992708632538, 'unit', 116, false, {'(pi^-4*sqrt(G))*(1+alpha*G^4)' '((pi^-4*sqrt(G))*(1+alpha*G^4))*(1-alpha*e^-4)' '(((pi^-4*sqrt(G))*(1+alpha*G^4))*(1-alpha*e^-4))*(1+alpha*e^-6)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'condenser_two_day', 0.0078904561464534125, 'unit', 116, false, {'(pi^-4*G^3)*(1+alpha*e^-4)' '((pi^-4*G^3)*(1+alpha*e^-4))*(1+alpha*e^-6)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'condenser_four_hour', 2.9554305405985954, 'positive', 116, false, {'(e+phi^-3)*(1+alpha*e^-3)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'condenser_twenty_hour', 1.8718105952496888, 'positive', 116, false, {'(G+sqrt(G))*(1-alpha*phi^-5)' '((G+sqrt(G))*(1-alpha*phi^-5))*(1+alpha*pi^-6)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'condenser_linear_scale', 0.66943295008216952, 'positive', 116, false, {'(phi^-5+gamma)*(1+alpha*phi^-2)' '((phi^-5+gamma)*(1+alpha*phi^-2))*(1+alpha*gamma^6)' '(((phi^-5+gamma)*(1+alpha*phi^-2))*(1+alpha*gamma^6))*(1+alpha*pi^-6)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'thermionic_delivered', 0.15151515151515152, 'positive', 116, false, {'(e*phi^-6)*(1+alpha*pi^-3)' '((e*phi^-6)*(1+alpha*pi^-3))*(1-alpha*pi^-5)' '(((e*phi^-6)*(1+alpha*pi^-3))*(1-alpha*pi^-5))*(1-alpha*pi^-6)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'thermionic_enrichment_return', 114.02750610678807, 'positive', 116, false, {'(e^5*G^3)*(1-alpha*pi^-3)' '((e^5*G^3)*(1-alpha*pi^-3))*(1+alpha*pi^-6)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'thermionic_headroom', 5.9999999999999991, 'positive', 116, false, {'(gamma^-2+gamma^-2)*(1-alpha*gamma^5)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'thermionic_ceiling_ratio', 0.90909090909090917, 'positive', 116, false, {'(gamma+gamma^2)*(1-alpha*gamma^3)' '((gamma+gamma^2)*(1-alpha*gamma^3))*(1-alpha*pi^-5)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'filed_margin', 0.0035714285714286143, 'unit', 116, false, {'(phi^-6*gamma^5)*(1+alpha*e^-4)' '((phi^-6*gamma^5)*(1+alpha*e^-4))*(1+alpha*e^-5)' '(((phi^-6*gamma^5)*(1+alpha*e^-4))*(1+alpha*e^-5))*(1-alpha*pi^-6)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'filed_rate', 50, 'positive', 116, false, {'(e^4*G)*(1-alpha*pi^-3)' '((e^4*G)*(1-alpha*pi^-3))*(1+alpha*pi^-5)' '(((e^4*G)*(1-alpha*pi^-3))*(1+alpha*pi^-5))*(1+alpha*pi^-6)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'filed_day_share', 0.0027397260273972603, 'unit', 116, false, {'(pi^-5*G^2)*(1-alpha*phi^-5)' '((pi^-5*G^2)*(1-alpha*phi^-5))*(1-alpha*e^-5)' '(((pi^-5*G^2)*(1-alpha*phi^-5))*(1-alpha*e^-5))*(1+alpha*pi^-6)' '((((pi^-5*G^2)*(1-alpha*phi^-5))*(1-alpha*e^-5))*(1+alpha*pi^-6))*(1+alpha*pi^-6)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'filed_cooldown_share', 0.0022831050228310501, 'unit', 116, false, {'(e^-6*G)*(1+alpha*sqrt(gamma))' '((e^-6*G)*(1+alpha*sqrt(gamma)))*(1+alpha*pi^-5)' '(((e^-6*G)*(1+alpha*sqrt(gamma)))*(1+alpha*pi^-5))*(1+alpha*pi^-6)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'filed_offset_ratio', 2000000, 'positive', 116, false, {'(pi^6*pi^6)*(1+alpha*e^5)' '((pi^6*pi^6)*(1+alpha*e^5))*(1+alpha*gamma^-3)' '(((pi^6*pi^6)*(1+alpha*e^5))*(1+alpha*gamma^-3))*(1+alpha*gamma^4)' '((((pi^6*pi^6)*(1+alpha*e^5))*(1+alpha*gamma^-3))*(1+alpha*gamma^4))*(1+alpha*pi^-5)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'filed_tube_ratio', 0.97331865859556077, 'positive', 116, false, {'(phi-G^5)*(1+alpha*e^-5)' '((phi-G^5)*(1+alpha*e^-5))*(1-alpha*pi^-6)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'filed_cooldown_ratio', 0.74063507354848279, 'positive', 116, false, {'(gamma^6+G^4)*(1-alpha*e^-3)' '((gamma^6+G^4)*(1-alpha*e^-3))*(1+alpha*e^-6)' '(((gamma^6+G^4)*(1-alpha*e^-3))*(1+alpha*e^-6))*(1-alpha*pi^-6)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'reference_kilogram', 0.98323377340494955, 'positive', 116, false, {'(phi^-4+G^2)*(1-alpha*phi^-3)' '((phi^-4+G^2)*(1-alpha*phi^-3))*(1+alpha*e^-5)' '(((phi^-4+G^2)*(1-alpha*phi^-3))*(1+alpha*e^-5))*(1-alpha*pi^-6)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'reference_core_ratio', 0.85371000000000008, 'positive', 116, false, {'(sqrt(pi)-G)*(1-alpha*phi^-2)' '((sqrt(pi)-G)*(1-alpha*phi^-2))*(1-alpha*gamma^5)' '(((sqrt(pi)-G)*(1-alpha*phi^-2))*(1-alpha*gamma^5))*(1+alpha*pi^-6)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'built_built_step', 1.0909090909090908, 'positive', 116, false, {'(G*G^-2)*(1-alpha*pi^-2)' '((G*G^-2)*(1-alpha*pi^-2))*(1-alpha*pi^-5)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'built_carnot_step', 1.2400793650793651, 'positive', 116, false, {'(e^-3+G^-2)*(1-alpha*gamma^3)' '((e^-3+G^-2)*(1-alpha*gamma^3))*(1+alpha*e^-4)' '(((e^-3+G^-2)*(1-alpha*gamma^3))*(1+alpha*e^-4))*(1-alpha*pi^-5)' '((((e^-3+G^-2)*(1-alpha*gamma^3))*(1+alpha*e^-4))*(1-alpha*pi^-5))*(1-alpha*pi^-6)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'tails_percent', 0.22, 'positive', 116, false, {'(phi^-2*gamma)*(1-alpha*gamma^2)' '((phi^-2*gamma)*(1-alpha*gamma^2))*(1+alpha*gamma^6)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'tails_fraction', 0.0022, 'unit', 116, false, {'(pi^-5-pi^-6)*(1-alpha*G^-6)' '((pi^-5-pi^-6)*(1-alpha*G^-6))*(1-alpha*pi^-5)' '(((pi^-5-pi^-6)*(1-alpha*G^-6))*(1-alpha*pi^-5))*(1-alpha*pi^-6)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'feed_multiple', 8.6831275720164616, 'positive', 116, false, {'(pi^2-G^-2)*(1+alpha*phi^-5)' '((pi^2-G^-2)*(1+alpha*phi^-5))*(1-alpha*pi^-5)' '(((pi^2-G^-2)*(1+alpha*phi^-5))*(1-alpha*pi^-5))*(1-alpha*pi^-6)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'feed_per_twh', 24.086757990867579, 'positive', 116, false, {'(pi^3-phi^4)*(1-alpha*phi^-2)' '((pi^3-phi^4)*(1-alpha*phi^-2))*(1+alpha*pi^-4)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'feed_concentrate', 0.84738955823293172, 'positive', 116, false, {'(e^-5+G^2)*(1+alpha*phi^-3)' '((e^-5+G^2)*(1+alpha*phi^-3))*(1+alpha*pi^-3)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'feed_hex', 0.67628205128205132, 'positive', 116, false, {'(pi^-3+G^5)*(1-alpha*phi^-4)' '((pi^-3+G^5)*(1-alpha*phi^-4))*(1-alpha*pi^-6)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'feed_enriched_hex', 0.67688022284122562, 'positive', 116, false, {'(pi^-3+G^5)*(1-alpha*pi^-3)' '((pi^-3+G^5)*(1-alpha*pi^-3))*(1+alpha*e^-5)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'feed_oxide', 0.88043478260869568, 'positive', 116, false, {'(sqrt(e)-G^3)*(1+alpha*pi^-3)' '((sqrt(e)-G^3)*(1+alpha*pi^-3))*(1-alpha*pi^-6)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'feed_core', 0.33287671232876714, 'positive', 116, false, {'(gamma*gamma)*(1-alpha*e^-2)' '((gamma*gamma)*(1-alpha*e^-2))*(1+alpha*pi^-4)' '(((gamma*gamma)*(1-alpha*e^-2))*(1+alpha*pi^-4))*(1+alpha*pi^-6)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'feed_transuranic', 0.010144927536231883, 'unit', 116, false, {'(pi*pi^-5)*(1-alpha*phi)' '((pi*pi^-5)*(1-alpha*phi))*(1+alpha*e^-6)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'ledger_quantum', 200, 'positive', 116, false, {'(e^2*gamma^-6)*(1+alpha*phi^-4)' '((e^2*gamma^-6)*(1+alpha*phi^-4))*(1+alpha*e^-6)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'aside_compare_concentrate', 0.84895833333333337, 'positive', 116, true, {'(pi^-4+G^2)*(1-alpha*e^-3)' '((pi^-4+G^2)*(1-alpha*e^-3))*(1+alpha*pi^-6)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'aside_compare_hex', 0.67634854771784236, 'positive', 116, true, {'(pi^-3+G^5)*(1-alpha*e^-2)' '((pi^-3+G^5)*(1-alpha*e^-2))*(1+alpha*e^-6)' '(((pi^-3+G^5)*(1-alpha*e^-2))*(1+alpha*e^-6))*(1-alpha*pi^-6)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'aside_compare_enriched_hex', 0.67400000000000004, 'positive', 116, true, {'(sqrt(G)*G^4)*(1+alpha*gamma^5)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'aside_compare_oxide', 0.88219895287958117, 'positive', 116, true, {'(phi^-3+G^5)*(1+alpha*gamma^3)' '((phi^-3+G^5)*(1+alpha*gamma^3))*(1+alpha*e^-4)' '(((phi^-3+G^5)*(1+alpha*gamma^3))*(1+alpha*e^-4))*(1+alpha*pi^-5)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'aside_compare_per_twh', 18.607305936073061, 'positive', 116, true, {'(gamma^-2+gamma^-5)*(1-alpha*e^-5)' '((gamma^-2+gamma^-5)*(1-alpha*e^-5))*(1+alpha*pi^-6)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'aside_compare_transuranic', 0.010471204188481676, 'unit', 116, true, {'(e^-5*G^-5)*(1+alpha*phi^-3)' '((e^-5*G^-5)*(1+alpha*phi^-3))*(1+alpha*gamma^6)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'aside_compare_used_oxide', 0.95811518324607325, 'positive', 116, true, {'(pi^-6+sqrt(G))*(1+alpha*e^-6)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'aside_red_tails', 0.25, 'positive', 116, false, {'(gamma^3*G^-3)*(1-alpha*e^-2)' '((gamma^3*G^-3)*(1-alpha*e^-2))*(1-alpha*e^-6)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'aside_red_demand', 175, 'positive', 116, false, {'(e^5+gamma^-6)*(1-alpha*gamma^2)' '((e^5+gamma^-6)*(1-alpha*gamma^2))*(1-alpha*e^-4)' '(((e^5+gamma^-6)*(1-alpha*gamma^2))*(1-alpha*e^-4))*(1-alpha*pi^-6)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'aside_red_lower_demand', 160, 'positive', 116, false, {'(e^5+phi^5)*(1+alpha*phi^-2)' '((e^5+phi^5)*(1+alpha*phi^-2))*(1+alpha*e^-3)' '(((e^5+phi^5)*(1+alpha*phi^-2))*(1+alpha*e^-3))*(1-alpha*e^-5)' '((((e^5+phi^5)*(1+alpha*phi^-2))*(1+alpha*e^-3))*(1-alpha*e^-5))*(1+alpha*pi^-6)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'aside_red_oxide', 206, 'positive', 116, false, {'(e^4+e^5)*(1+alpha*sqrt(pi))' '((e^4+e^5)*(1+alpha*sqrt(pi)))*(1+alpha*phi^-3)' '(((e^4+e^5)*(1+alpha*sqrt(pi)))*(1+alpha*phi^-3))*(1+alpha*e^-5)' '((((e^4+e^5)*(1+alpha*sqrt(pi)))*(1+alpha*phi^-3))*(1+alpha*e^-5))*(1-alpha*pi^-6)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'aside_red_lower_oxide', 189, 'positive', 116, false, {'(pi^5/phi)*(1-alpha*phi^-5)' '((pi^5/phi)*(1-alpha*phi^-5))*(1-alpha*pi^-5)' '(((pi^5/phi)*(1-alpha*phi^-5))*(1-alpha*pi^-5))*(1-alpha*pi^-6)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'aside_candu_uranium', 98.799999999999997, 'positive', 116, false, {'(pi^4+G^-4)*(1-alpha*gamma^6)' '((pi^4+G^-4)*(1-alpha*gamma^6))*(1-alpha*pi^-5)' '(((pi^4+G^-4)*(1-alpha*gamma^6))*(1-alpha*pi^-5))*(1-alpha*pi^-6)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'aside_candu_plutonium', 0.38, 'positive', 116, false, {'(G^5*G^6)*(1-alpha*phi^-3)' '((G^5*G^6)*(1-alpha*phi^-3))*(1-alpha*e^-3)' '(((G^5*G^6)*(1-alpha*phi^-3))*(1-alpha*e^-3))*(1+alpha*e^-5)' '((((G^5*G^6)*(1-alpha*phi^-3))*(1-alpha*e^-3))*(1+alpha*e^-5))*(1-alpha*pi^-6)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'aside_candu_efficiency', 0.31, 'positive', 116, false, {'(phi^-3/sqrt(gamma))*(1-alpha*gamma^2)' '((phi^-3/sqrt(gamma))*(1-alpha*gamma^2))*(1+alpha*e^-4)' '(((phi^-3/sqrt(gamma))*(1-alpha*gamma^2))*(1+alpha*e^-4))*(1-alpha*e^-6)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'aside_candu_demand', 157, 'positive', 116, false, {'(e^5+gamma^-4)*(1-alpha*phi^-2)' '((e^5+gamma^-4)*(1-alpha*phi^-2))*(1+alpha*e^-4)' '(((e^5+gamma^-4)*(1-alpha*phi^-2))*(1+alpha*e^-4))*(1-alpha*pi^-5)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'aside_candu_capacity', 142, 'positive', 116, false, {'(e^5*sqrt(G))*(1-alpha*gamma^6)' '((e^5*sqrt(G))*(1-alpha*gamma^6))*(1-alpha*e^-6)'});
fails = check_close(fid, fails, S, base_names, base_values, helium, 'year_share', 0.00045662100456621003, 'unit', 116, false, {'(pi^-5*e^-2)*(1+alpha*phi^3)' '((pi^-5*e^-2)*(1+alpha*phi^3))*(1+alpha*gamma^3)' '(((pi^-5*e^-2)*(1+alpha*phi^3))*(1+alpha*gamma^3))*(1+alpha*e^-4)' '((((pi^-5*e^-2)*(1+alpha*phi^3))*(1+alpha*gamma^3))*(1+alpha*e^-4))*(1+alpha*pi^-6)' '(((((pi^-5*e^-2)*(1+alpha*phi^3))*(1+alpha*gamma^3))*(1+alpha*e^-4))*(1+alpha*pi^-6))*(1+alpha*pi^-6)'});
fails = check_pair(fid, fails, S, base_names, base_values, helium, 'removal_one_day_pair', 1.7697019567375738, 'positive', '(sqrt(pi)-e^-6)*(1-alpha*e^-4)', '((sqrt(pi)-e^-6)*(1-alpha*e^-4))*(1-alpha*e^-6)', '((sqrt(pi)-e^-6)*(1-alpha*e^-4))*(1-alpha*pi^-5)');
fails = check_pair(fid, fails, S, base_names, base_values, helium, 'margin_proportion_pair', 0.0027000000000000001, 'unit', '((e^-6/G)*(1-alpha*gamma^2))*(1+alpha*e^-4)', '(((e^-6/G)*(1-alpha*gamma^2))*(1+alpha*e^-4))*(1+alpha*e^-6)', '(((e^-6/G)*(1-alpha*gamma^2))*(1+alpha*e^-4))*(1+alpha*pi^-5)');
fails = check_pair(fid, fails, S, base_names, base_values, helium, 'stack_fraction_pair', 0.044999999999999998, 'unit', '((e^-2*gamma^2)*(1-alpha*phi^-3))*(1-alpha*gamma^6)', '(((e^-2*gamma^2)*(1-alpha*phi^-3))*(1-alpha*gamma^6))*(1-alpha*e^-6)', '(((e^-2*gamma^2)*(1-alpha*phi^-3))*(1-alpha*gamma^6))*(1-alpha*pi^-5)');
fails = check_pair(fid, fails, S, base_names, base_values, helium, 'filed_rejected', 0.5535714285714286, 'positive', '((gamma*sqrt(G))*(1+alpha*phi^-3))*(1+alpha*e^-3)', '(((gamma*sqrt(G))*(1+alpha*phi^-3))*(1+alpha*e^-3))*(1-alpha*e^-6)', '(((gamma*sqrt(G))*(1+alpha*phi^-3))*(1+alpha*e^-3))*(1-alpha*pi^-5)');
fprintf(fid, 'closes %d\n', 97);
fprintf(fid, 'pairs %d\n', 4);
fprintf(fid, 'miss 0\n');
fprintf(fid, 'fails %d\n', fails);
fclose(fid);
if fails > 0
  error('fission stamp failed: %d', fails);
end
endfunction


function helium = helium_match(S)
bare = S.pi / S.gamma ^ 4;
central = 28.2956624;
quotient = (bare - central) / bare / S.alpha ^ 2;
witness = S.pi + S.gamma / S.e;
helium = abs((witness - quotient) / quotient);
endfunction

function [names, values] = base_pieces(S)
atoms = {'pi', 'e', 'phi', 'gamma', 'G'};
avals = [S.pi, S.e, S.phi, S.gamma, S.G];
names = cell(1, 60);
values = zeros(1, 60);
cursor = 1;
for index = 1:5
  names{cursor} = atoms{index};
  values(cursor) = avals(index);
  cursor = cursor + 1;
  names{cursor} = ['sqrt(' atoms{index} ')'];
  values(cursor) = sqrt(avals(index));
  cursor = cursor + 1;
  for power = 2:6
    names{cursor} = sprintf('%s^%d', atoms{index}, power);
    values(cursor) = avals(index) ^ power;
    cursor = cursor + 1;
    names{cursor} = sprintf('%s^-%d', atoms{index}, power);
    values(cursor) = avals(index) ^ (-power);
    cursor = cursor + 1;
  end
end
endfunction

function [names, values] = dressings(anchor_name, anchor, alpha, base_names, base_values, mode)
names = {};
values = [];
for index = 1:length(base_names)
  for sign = [1, -1]
    number = anchor * (1 + sign * alpha * base_values(index));
    if ~isfinite(number)
      continue
    end
    if strcmp(mode, 'unit')
      keep = number > 0 && number < 1;
    else
      keep = number > 0;
    end
    if keep
      if sign < 0
        mark = '-';
      else
        mark = '+';
      end
      names{end + 1} = ['(' anchor_name ')*(1' mark 'alpha*' base_names{index} ')'];
      values(end + 1) = number;
    end
  end
end
endfunction

function spacing = median_spacing(values)
rounded = zeros(size(values));
for index = 1:length(values)
  rounded(index) = py_round12(values(index));
end
vals = unique(rounded);
vals = sort(vals);
spacings = [];
for index = 1:(length(vals) - 1)
  left = vals(index);
  right = vals(index + 1);
  if right > left
    spacings(end + 1) = (right - left) / right;
  end
end
spacings = sort(spacings);
spacing = spacings(floor(length(spacings) / 2) + 1);
endfunction

function value = py_round12(number)
scale = 1e12;
scaled = number * scale;
whole = floor(scaled);
frac = scaled - whole;
if frac > 0.5 || (frac == 0.5 && mod(whole, 2) ~= 0)
  whole = whole + 1;
end
value = whole / scale;
endfunction

function value = evaluate(text, S)
[value, cursor] = parse_expr(text, 1, S);
if cursor ~= length(text) + 1
  error('unparsed expression');
end
endfunction

function [value, cursor] = parse_expr(text, cursor, S)
[value, cursor] = parse_term(text, cursor, S);
while cursor <= length(text) && (text(cursor) == '+' || text(cursor) == '-')
  op = text(cursor);
  [right, cursor] = parse_term(text, cursor + 1, S);
  if op == '+'
    value = value + right;
  else
    value = value - right;
  end
end
endfunction

function [value, cursor] = parse_term(text, cursor, S)
[value, cursor] = parse_power(text, cursor, S);
while cursor <= length(text) && (text(cursor) == '*' || text(cursor) == '/')
  op = text(cursor);
  [right, cursor] = parse_power(text, cursor + 1, S);
  if op == '*'
    value = value * right;
  else
    value = value / right;
  end
end
endfunction

function [value, cursor] = parse_power(text, cursor, S)
[value, cursor] = parse_unary(text, cursor, S);
if cursor <= length(text) && text(cursor) == '^'
  [exponent, cursor] = parse_unary(text, cursor + 1, S);
  value = value ^ exponent;
end
endfunction

function [value, cursor] = parse_unary(text, cursor, S)
if cursor <= length(text) && text(cursor) == '-'
  [value, cursor] = parse_unary(text, cursor + 1, S);
  value = -value;
  return
end
[value, cursor] = parse_primary(text, cursor, S);
endfunction

function [value, cursor] = parse_primary(text, cursor, S)
if text(cursor) == '('
  [value, cursor] = parse_expr(text, cursor + 1, S);
  if text(cursor) ~= ')'
    error('missing parenthesis');
  end
  cursor = cursor + 1;
  return
end
if starts_token(text, cursor, 'sqrt(')
  [value, cursor] = parse_expr(text, cursor + 5, S);
  cursor = cursor + 1;
  value = sqrt(value);
  return
end
if text(cursor) >= '0' && text(cursor) <= '9'
  start = cursor;
  while cursor <= length(text) && text(cursor) >= '0' && text(cursor) <= '9'
    cursor = cursor + 1;
  end
  value = str2double(text(start:cursor - 1));
  return
end
for name = {'alpha', 'gamma', 'phi', 'pi', 'e', 'G'}
  token = name{1};
  if starts_token(text, cursor, token)
    cursor = cursor + length(token);
    value = S.(token);
    return
  end
end
error('bad token');
endfunction

function yes = starts_token(text, cursor, token)
last = cursor + length(token) - 1;
yes = last <= length(text) && strcmp(text(cursor:last), token);
if ~yes
  return
end
if strcmp(token, 'sqrt(')
  return
end
if last == length(text)
  return
end
next = text(last + 1);
yes = next == '+' || next == '-' || next == '*' || next == '/' || next == '^' || next == ')';
endfunction

function anchor = parent_anchor(expr)
cut = 0;
tokens = {'*(1+alpha*', '*(1-alpha*'};
for index = 1:2
  hits = strfind(expr, tokens{index});
  if ~isempty(hits)
    cut = max(cut, hits(end) - 1);
  end
end
anchor = expr(2:cut - 1);
endfunction

function fails = check_close(fid, fails, S, base_names, base_values, helium, ident, target, mode, count, unrun, steps)
anchor = parent_anchor(steps{1});
for index = 1:length(steps)
  [names, values] = dressings(anchor, evaluate(anchor, S), S.alpha, base_names, base_values, mode);
  gaps = abs(values - target) ./ abs(target);
  [~, order] = sortrows([gaps(:), (1:length(gaps))']);
  best = order(1);
  inside = sum(gaps <= helium);
  spacing = median_spacing(values);
  expected_inside = 0;
  if index == length(steps)
    expected_inside = 1;
  end
  ok = strcmp(names{best}, steps{index}) && inside == expected_inside && length(names) == count;
  ok = ok && length(unique(arrayfun(@py_round12, values))) == count;
  if count == 116
    ok = ok && abs(spacing - 0.000511364) < 2e-9 && spacing > helium;
  else
    ok = ok && count == 110 && abs(spacing - 0.00043369787200611595) < 2e-9 && spacing > helium;
  end
  if ~ok
    fails = fails + 1;
    fprintf(fid, 'FAIL %s step %d inside %d count %d spacing %.12e nearest %s\n', ident, index, inside, length(names), spacing, names{best});
  end
  anchor = steps{index};
end
closing_gap = abs(evaluate(steps{end}, S) - target) / abs(target);
if closing_gap > helium
  fails = fails + 1;
  fprintf(fid, 'FAIL %s closing gap %.12e\n', ident, closing_gap);
end
status = 'named';
if unrun
  status = 'unrun';
end
fprintf(fid, '%s %s gap %.6e factors %d count %d\n', ident, status, closing_gap, length(steps), count);
endfunction

function fails = check_pair(fid, fails, S, base_names, base_values, helium, ident, target, mode, parent, e_name, pi_name)
parent_gap = abs(evaluate(parent, S) - target) / abs(target);
e_gap = abs(evaluate(e_name, S) - target) / abs(target);
pi_gap = abs(evaluate(pi_name, S) - target) / abs(target);
[names, values] = dressings(parent, evaluate(parent, S), S.alpha, base_names, base_values, mode);
gaps = abs(values - target) ./ abs(target);
inside_names = names(gaps <= helium);
spacing = median_spacing(values);
ok = parent_gap > helium && e_gap <= helium && pi_gap <= helium;
ok = ok && length(inside_names) == 2;
ok = ok && any(strcmp(inside_names, e_name)) && any(strcmp(inside_names, pi_name));
ok = ok && length(names) == 116 && abs(spacing - 0.000511364) < 2e-9;
if ~ok
  fails = fails + 1;
  fprintf(fid, 'FAIL %s pair inside %d parent %.6e e %.6e pi %.6e\n', ident, length(inside_names), parent_gap, e_gap, pi_gap);
end
for child = {e_name, pi_name}
  child_name = child{1};
  child_gap = abs(evaluate(child_name, S) - target) / abs(target);
  [next_names, next_values] = dressings(child_name, evaluate(child_name, S), S.alpha, base_names, base_values, mode);
  next_gaps = abs(next_values - target) ./ abs(target);
  if min(next_gaps) <= child_gap || length(next_names) ~= 116
    fails = fails + 1;
    fprintf(fid, 'FAIL %s next factor tightened\n', ident);
  end
end
fprintf(fid, '%s held e %.6e pi %.6e\n', ident, e_gap, pi_gap);
endfunction


