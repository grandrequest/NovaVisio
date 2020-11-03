cellTypes = {'ON S' 'OFF T' 'OFF S'};
cellLabels = containers.Map;

% 20200716
cellInfo = containers.Map;
cellInfo('Cell11') = {'ON S' 'KO'};
cellInfo('Cell2') = {'ON S' 'WT'};
cellInfo('Cell4') = {'ON S' 'WT'};
cellInfo('Cell5') = {'ON S' 'WT'};
cellInfo('Cell7') = {'ON S' 'KO'};
cellInfo('Cell3') = {'OFF T' 'WT'};
cellInfo('Cell6') = {'OFF T' 'WT'};
cellInfo('Cell8')  = {'OFF T' 'KO'};
cellInfo('Cell12') = {'OFF T' 'KO'};
cellInfo('Cell1')= {'OFF S' 'WT'};
cellInfo('Cell9')  = {'OFF S' 'KO'};
cellInfo('Cell10')  = {'OFF S' 'KO'};

cellLabels('20200716') = cellInfo;

% 20200728
cellInfo = containers.Map;
cellInfo('Cell2') = {'OFF S' 'WT'};
cellInfo('Cell4') = {'ON S' 'WT'};
cellInfo('Cell5') = {'OFF S' 'WT'};
cellInfo('Cell6') = {'ON S' 'WT'};
cellInfo('Cell1') = {'ON S' 'WT'};

cellLabels('20200728') = cellInfo;


%20200721

cellInfo = containers.Map;
cellInfo('Cell1')= {"ON S" "KO"};
cellInfo('Cell2')= {"ON S" "KO"};
cellInfo('Cell4')= {"ON S" "KO"};
cellInfo('Cell7')= {"ON S" "KO"};
cellInfo('Cell8')= {"OFF T" "KO"};
cellInfo('Cell12')= {"OFF T" "KO"};
cellInfo('Cell3')= {"OFF S" "KO"};
cellInfo('Cell9')= {"OFF S" "KO"};
cellInfo('Cell10')= {"OFF S" "KO"};

cellLabels('20200721') = cellInfo;

%20201022

cellInfo = containers.Map;
cellInfo('Cell1')= {"ON S" "KO"};
cellInfo('Cell2')= {"ON S" "KO"};

cellLabels('20201022') = cellInfo;

%20201007

cellInfo = containers.Map;
cellInfo('Cell1')= {"ON S" "KO"};
cellInfo('Cell2')= {"ON S" "KO"};

cellLabels('20201007') = cellInfo;

%20200924

cellInfo = containers.Map;
cellInfo('Cell1')= {"ON S" "KO"};
cellInfo('Cell2')= {"ON S" "KO"};
cellInfo('Cell3')= {"ON S" "KO"};
cellInfo('Cell4')= {"ON S" "WT"};
cellInfo('Cell5')= {"ON S" "WT"};
cellInfo('Cell6')= {"OFF S" "WT"};
cellInfo('Cell7')= {"ON S" "WT"};

cellLabels('20200924') = cellInfo;

%20200723

cellInfo = containers.Map;
cellInfo('Cell1')= {"OFF S" "WT"};
cellInfo('Cell2')= {"ON S" "WT"};
cellInfo('Cell3')= {"OFF S" "WT"};
cellInfo('Cell5')= {"ON S" "WT"};
cellInfo('Cell6')= {"OFF S" "WT"};
cellInfo('Cell7')= {"ON S" "WT"};
cellInfo('Cell8')= {"OFF T" "WT"};

cellLabels('20200723') = cellInfo;