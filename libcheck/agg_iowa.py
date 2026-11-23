import sys, pandas as pd
parts=[]
for i,ch in enumerate(pd.read_csv(sys.stdin, usecols=['Date','County Number','Sale (Dollars)'], dtype={'Date':str,'County Number':'float64','Sale (Dollars)':'float64'}, chunksize=2_000_000)):
    g=ch.groupby(['Date','County Number'])['Sale (Dollars)'].sum()
    parts.append(g)
    if i%5==0: print(i, len(ch), flush=True)
s=pd.concat(parts).groupby(level=[0,1]).sum().reset_index()
s.columns=['date','county','sales']
s['date']=pd.to_datetime(s['date'],format='%m/%d/%Y')
s.to_parquet('iowa_daily_county.parquet')
print('done', s.shape, s.date.min(), s.date.max(), s.county.nunique())
