%initCobraToolbox;
TAC125Model = readCbModel('iMF721_SMBL3FBC2.xml');
%TAC125Model = checkCobraModelUnique(TAC125Model, 'T')
FBAsolutionAsIs = optimizeCbModel(TAC125Model,'max');
GrowthAllEXOpen = FBAsolutionAsIs.f
EX_reactions= TAC125Model.rxns(~cellfun(@isempty, regexp(TAC125Model.rxns,'^EX_')));
NumberOfEXReactions = length(EX_reactions)
TAC125Model = changeRxnBounds(TAC125Model, EX_reactions, 0, 'l');
FBAsolutionNoEX = optimizeCbModel(TAC125Model,'max');
FBAsolutionAllClosed = FBAsolutionNoEX.f
RxnListExchangeSchatz = { 'EX_cpd00254_e' , 'EX_cpd00012_e' , 'EX_cpd00971_e' , 'EX_cpd00067_e' , 'EX_cpd00063_e' , 'EX_cpd00048_e' , 'EX_cpd00205_e' , 'EX_cpd00099_e' , 'EX_cpd00009_e' , 'EX_cpd00013_e' , 'EX_cpd10515_e',  'EX_cpd00209_e' };
TAC125Model = changeRxnBounds(TAC125Model, RxnListExchangeSchatz, -1000, 'l');

SchatzCompunds = strrep(RxnListExchangeSchatz, 'EX_', '');
SchatzCompunds = strrep(SchatzCompunds, '(e)','[e]');
SchatzCompunds = strrep(SchatzCompunds, '_c','c]');

SchatzMatched = TAC125Model.mets(ismember(TAC125Model.mets, SchatzCompunds'));
SchatzMatchedNames = TAC125Model.metNames(ismember(TAC125Model.mets, SchatzCompunds'));
SchatzMatchedLB = TAC125Model.lb(findRxnIDs(TAC125Model, RxnListExchangeSchatz));
SchatzMatchedlUB = TAC125Model.ub(findRxnIDs(TAC125Model, RxnListExchangeSchatz));

table(SchatzMatched, SchatzMatchedNames)
table(RxnListExchangeSchatz', SchatzMatchedlUB, SchatzMatchedLB)


LeuGrowthRateExp = 0.10;
GluGrowthRateExp = 0.34;
AlaGrowthRateExp = 0.25;
AspGrowthRateExp = 0.20;
GGGrowthRateExp = 1.1;


LeuUtilizationRate = 0.1/(5.4/5)/131.17*1000;
GluUtilizationRate = 0.34/(4.4/6.5)/147.13*1000;
AspUtilizationRate = 0.25/(3.6/5.1)/133.11*1000;
AlaUtilizationRate = 0.20/(2.1/3.4)/89.09*1000;

TAC125Model = changeRxnBounds(TAC125Model,  'EX_cpd00107_e', -LeuUtilizationRate, 'l');
FBAsolutionSchatzLeu = optimizeCbModel(TAC125Model,'max');
SimulationLeu = FBAsolutionSchatzLeu.f;  
TAC125Model = changeRxnBounds(TAC125Model,  'EX_cpd00107_e', 0, 'l');

TAC125Model = changeRxnBounds(TAC125Model,  'EX_cpd00023_e', -3.3, 'l');
FBAsolutionSchatzGlu = optimizeCbModel(TAC125Model,'max');
SimulationGlu = FBAsolutionSchatzGlu.f;  
TAC125Model = changeRxnBounds(TAC125Model,  'EX_cpd00023_e', 0, 'l');

TAC125Model = changeRxnBounds(TAC125Model,  'EX_cpd00035_e', -AlaUtilizationRate, 'l');
FBAsolutionSchatzAla = optimizeCbModel(TAC125Model,'max');
SimulationAla = FBAsolutionSchatzAla.f;  
TAC125Model = changeRxnBounds(TAC125Model,  'EX_cpd00035_e', 0, 'l');

TAC125Model = changeRxnBounds(TAC125Model,  'EX_cpd00041_e', -AspUtilizationRate , 'l');
FBAsolutionSchatzAsp = optimizeCbModel(TAC125Model,'max');
SimulationAsp = FBAsolutionSchatzAsp.f  
TAC125Model = changeRxnBounds(TAC125Model,  'EX_cpd00041_e', 0, 'l');


table(LeuUtilizationRate, SimulationLeu, LeuGrowthRateExp)
table(GluUtilizationRate, SimulationGlu, GluGrowthRateExp)
table(AlaUtilizationRate, SimulationAla, AlaGrowthRateExp)
table(AspUtilizationRate, SimulationAsp, AspGrowthRateExp)

%Set GG simulation

TAC125Model = changeRxnBounds(TAC125Model,  'EX_cpd00023_e', -.1, 'l');
TAC125Model = changeRxnBounds(TAC125Model,  'EX_cpd00222_e', -.1, 'l');
FBAsolutionSchatzGG= optimizeCbModel(TAC125Model,'max');
SimulationGG = FBAsolutionSchatzGG.f
ErrorsInExp = [0.01, 0.03, 0.02, 0.25, 0.03]
ModelGR = [SimulationLeu, SimulationGlu, SimulationAla, SimulationAsp, SimulationGG ]
ExperimentalGR = [LeuGrowthRateExp, GluGrowthRateExp, AlaGrowthRateExp, AspGrowthRateExp, GGGrowthRateExp ]
Combined = [ModelGR(:), ExperimentalGR(:)]
bar(Combined)

ylabel('Growth rate (h^-^1)')
xticklabels({'Leu', 'Glu', 'Ala', 'Asp', 'GG'})
set(gca,'FontName','Times','fontsize',24)
legend('Experimental','Model', 'Location','northwest')

AAExchange = { 'EX_cpd00254_e' , 'EX_cpd00012_e' , 'EX_cpd00971_e' , 'EX_cpd00067_e' , 'EX_cpd00063_e' , 'EX_cpd00048_e' , 'EX_cpd00205_e' , 'EX_cpd00099_e' , 'EX_cpd00009_e' , 'EX_cpd00013_e' , 'EX_cpd10515_e',  'EX_cpd00209_e' };

RxnListExchangeAA = {'EX_cpd00119_e' , 'EX_cpd00322_e' , 'EX_cpd00107_e' , 'EX_cpd00039_e' ,'EX_cpd00060_e' , 'EX_cpd00066_e' ,  'EX_cpd00161_e' , 'EX_cpd00065_e' , 'EX_cpd00156_e', 'EX_cpd00035_e' , 'EX_cpd00051_e' , 'EX_cpd00132_e'  ,  'EX_cpd00041_e' , 'EX_cpd00084_e', 'EX_cpd00023_e' , 'EX_cpd00053_e' , 'EX_cpd00033_e' , 'EX_cpd00064_e' , 'EX_cpd00129_e'  'EX_cpd00054_e' , 'EX_cpd00069_e'}

for i=1:length(RxnListExchangeAA)

    TAC125ModelAA = changeRxnBounds(TAC125Model,  RxnListExchangeAA(i), -1, 'l');
    optimizeCbModel(TAC125ModelAA)
    
end