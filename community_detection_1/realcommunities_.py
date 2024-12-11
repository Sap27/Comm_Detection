# Function below is specifically written for LFR benchmarks that come with a community.dat file (not yet tried communities with overlapping nodes )
import pandas as pd
def groundtruth(community_path):
    community=open(community_path)
    nodeid=[]
    communityid=[]
    for i in community:
        l=len(i.split('\t')[1].split(' '))-1
        for j in range(l):
            nodeid.append(int(i.split('\t')[0]))
            communityid.append(int(i.split('\t')[1].split(' ')[j]))
    dict_={'nodeid':nodeid,'communityid':communityid}
    df=pd.DataFrame(dict_)
    grouped=df.groupby(df.communityid)
    realcommunities=[]
    for i in range(max(communityid)):
        df_new=grouped.get_group(i+1)
        realcommunities.append(list(df_new['nodeid']))
    return realcommunities