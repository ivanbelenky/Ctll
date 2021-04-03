import matplotlib.pyplot as plt

import overpy
api = overpy.Overpass()

result=api.query("area['name:en'='Denmark']->.country;rel['name:en'='Denmark']['type'='boundary']['admin_level'='2'];(way(r)['maritime' != 'yes'](40,-10,70,80);way(area.country)['natural'='coastline'](40,-10,70,80););out geom;")

x=[]
y=[]
i=0
for way in result.ways:
    print(f"way {i} of {len(result.ways)}")
    if 'natural' in way.tags and way.tags['natural']=='coastline' and len(way.get_nodes(True))>0: #just a test
        i=i+1
        j=0
        for node in way.get_nodes(True):
            j+=1
            if j%100 != 0:
            	continue
            print (f'lon: {float(node.lon):3.4f}; lat: {float(node.lat):3.4f}')
            x.append(float(node.lon))
            y.append(float(node.lat))
plt.plot(x, y, 'o',label=str(way.id))
plt.show()