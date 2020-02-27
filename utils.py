"""
Utility functions to create list and network of possible supervisors and/or collaborators.
Author: Lorenzo Fabbri
"""

from Bio import Entrez
from pprint import pprint
import pandas as pd
from geotext import GeoText
from geopy.geocoders import Nominatim
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import folium

class Search:

    def __init__(self, email, retmax):
        self.query = ""
        self.countries = []
        self.email = email
        self.retmax = retmax
    
    def esearch(self, query):
        self.query = query
        Entrez.email = self.email

        handle = Entrez.esearch(db='pubmed', sort='relevance', 
                                retmax=self.retmax, retmode='xml', 
                                term=self.query)
        results = Entrez.read(handle)

        return results
    
    def efetch(self, ids_list):
        ids = ','.join(ids_list)
        Entrez.email = self.email

        handle = Entrez.efetch(db='pubmed', retmode='xml', 
                                id=ids)
        results = Entrez.read(handle)

        return results
    
    def search(self, query):
        results = self.esearch(query)
        ids_list = results['IdList']
        papers = self.efetch(ids_list)

        df_authors_papers = []
        for idx, paper in enumerate(papers['PubmedArticle']):
            try:
                article = paper['MedlineCitation']['Article']
                title = article['ArticleTitle']
                abstract = article['Abstract']
                authors = article['AuthorList']
            except:
                continue

            if authors:
                temp = []
                for author in authors:
                    try:
                        affiliation = author['AffiliationInfo'][0]['Affiliation']
                    except:
                        continue
                    name = author['ForeName']
                    surname = author['LastName']

                    temp_author = {'name': name, 
                                    'surname': surname, 
                                    'affiliation': affiliation}
                    temp.append(temp_author)
            
            df_authors_papers.append(pd.DataFrame(temp))
        
        df = pd.concat(df_authors_papers, ignore_index=True, 
                        sort=True)
        df.sort_values(by=['surname'], inplace=True)
        df.drop_duplicates(subset=['name', 'surname'], inplace=True)
        df.reset_index(inplace=True, drop=True)

        return df
    
    def extract_geo(self, df):
        #df.drop_duplicates(subset=['affiliation'], inplace=True)
        df.reset_index(inplace=True, drop=True)

        data = []
        for idx in df.index:
            affiliation = df.loc[df.index[idx], 'affiliation']
            city = GeoText(affiliation).cities
            if city:
                data.append([' '.join(city), 
                            affiliation, 
                            df.loc[df.index[idx], 'name'], 
                            df.loc[df.index[idx], 'surname']])
        
        return data
    
    def map(self, data):
        geolocator = Nominatim(timeout=10)

        map = Basemap(width=10000000, height=6000000, 
                        projection='cyl', resolution='i')
        map.drawcountries()
        map.drawcoastlines()
        map.drawmapboundary()
        map.fillcontinents(color='grey')

        for datum in data:
            city = datum[0]
            loc = geolocator.geocode(city)
            if not loc:
                print(f'Cannot locate {city}.')
                continue
            
            x, y = map(loc.longitude, loc.latitude)
            map.plot(x, y, marker='o', color='Red')
            plt.annotate(datum[2] + ' ' + datum[3], xy=(x, y))
        
        plt.show()
    
    def interactive_map(self, data):

        folium_map = folium.Map()
        geolocator = Nominatim(timeout=10)

        for datum in data:
            city = datum[0]
            loc = geolocator.geocode(city)
            if not loc:
                print(f'Cannot locate {city}.')
                continue
            
            lat = loc.latitude
            lon = loc.longitude
            txt = f'Affiliation: {datum[1]}<br> Name: {datum[2]}<br> Surname: {datum[3]}'
            folium.CircleMarker(location=(lat, lon), 
                                popup=txt, 
                                fill=True).add_to(folium_map)
        
        return folium_map
