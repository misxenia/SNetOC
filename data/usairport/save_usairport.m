clear 

%% http://toreopsahl.com/datasets/#usairports

fid = fopen('USairport_2010.txt');
A = textscan(fid, '%f%f%f');
fclose(fid);
G = sparse(A{1}, A{2}, A{3});

% airport codes
fid = fopen('USairport_2010_codes.txt');
A = textscan(fid, '%d%q');
code = A{2};
fclose(fid);

% %% Load BTS Transtats data
% % Downloaded from http://www.transtats.bts.gov/DL_SelectFields.asp?Table_ID=292
% % Filters: Geography=all; Year=2010; Months=all
% unzip('577903858_T_T100_MARKET_ALL_CARRIER.zip')
% fid = fopen('577903858_T_T100_MARKET_ALL_CARRIER_2010_All.csv');
% colnames = textscan(fid, repmat('%q',1,41), 1, 'delimiter', ',');
% colnames = [colnames{:}];
% %                  [flight][carrier         ][origin              ][destination         ][time  ]
% A = textscan(fid, '%f%f%f%f%q%f%q%q%q%q%q%f%f%f%f%f%q%q%q%q%q%q%q%f%f%f%f%q%q%q%q%q%q%q%f%f%f%f%f%q%q', 'delimiter', ',');
% fclose(fid);
% A = cell2struct(A, colnames, 2);
% 
% ids = unique([A.ORIGIN_AIRPORT_ID; A.DEST_AIRPORT_ID]);
% K = numel(ids);
% 
% [~, ind_org] = ismember(A.ORIGIN_AIRPORT_ID, ids);
% [~, ind_dest] = ismember(A.DEST_AIRPORT_ID, ids);
% G = logical(sparse(ind_org, ind_dest, A.PASSENGERS, K, K));

%% airports meta data
% Load airport information (incl. geolocations)
% Downloaded from http://www.transtats.bts.gov/DL_SelectFields.asp?Table_ID=288 (select all columns)

unzip('253021595_T_MASTER_CORD.zip')
fid = fopen('253021595_T_MASTER_CORD_All_All.csv', 'r');
colnames = textscan(fid, repmat('%q', 1,32), 1, 'Delimiter', ',');
colnames = [colnames{:}];
A = textscan(fid, '%f%f%q%q%q%f%f%q%q%q%q%q%f%f%q%f%f%f%q%f%f%f%f%q%f%f%f%q%q%q%f%f', 'Delimiter', ',', 'EndOfLine', '\r\n');
fclose(fid);
A = cell2struct(A, colnames, 2);
A.AIRPORT_START_DATE = datetime(A.AIRPORT_START_DATE, 'InputFormat', 'yyy-MM-dd');
A.AIRPORT_THRU_DATE = datetime(A.AIRPORT_THRU_DATE, 'InputFormat', 'yyy-MM-dd');
% For code requiring serial dates (datenum) instead of datetime, uncomment
% the following line(s) below to return the imported dates as datenum(s).
%A.AIRPORT_START_DATE = datenum(AIRPORT_START_DATE);
%A.AIRPORT_THRU_DATE = datenum(AIRPORT_THRU_DATE);

%[~, ind] = ismember(ids, A.AIRPORT_ID);


[~, ind] = ismember(code, A.AIRPORT);

meta.code = A.AIRPORT(ind);
meta.name = A.DISPLAY_AIRPORT_NAME(ind);
meta.city = A.DISPLAY_AIRPORT_CITY_NAME_FULL(ind);
meta.country = A.AIRPORT_COUNTRY_NAME(ind);
meta.country_code = A.AIRPORT_COUNTRY_CODE_ISO(ind);
meta.state = A.AIRPORT_STATE_NAME(ind);
meta.state_code = A.AIRPORT_STATE_CODE(ind);
meta.lat = A.LATITUDE(ind);
meta.lon = A.LONGITUDE(ind);

%%
G = G | G'; % make undirected graph
G = G-diag(diag(G)); % remove self loops (#120)

% remove nodes without edges (#281)
ok = sum(G)>0;
G = G(ok,ok);

fn = fieldnames(meta);
for i=1:numel(fn)
    meta.(fn{i}) = meta.(fn{i})(ok);
end
 
%% save data
save usairport.mat G meta
