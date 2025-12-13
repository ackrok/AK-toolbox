load('C:\Users\Anya\Desktop\FP_LOCAL\rev_C\achda_beh_cannulaSNc_saline+muscimol');

compare = cannula;
compare(1).lbl = compare(1).inf; compare(2).lbl = compare(2).inf;
winInf = [20 40]; winInf = winInf.*(60*50);
for z = 1:length(compare)
    for x = 1:length(compare(z).s)
        rmv = find(compare(z).s(x).onRest < winInf(1));
        rmv = [rmv; find(compare(z).s(x).onRest > winInf(2))];
        compare(z).s(x).onRest(rmv) = [];
        compare(z).s(x).offRest(rmv) = [];
    end
end

%%
plot_revD_compareImmTrPk
subplot(2,3,1); ylim([0 1]);
subplot(2,3,4); ylim([0 1]);
