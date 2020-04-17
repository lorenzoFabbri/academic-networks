"""
Utility classes.
Author: Lorenzo Fabbri
"""

import pandas as pd
from Bio import Entrez
from geotext import GeoText
from geopy.geocoders import Nominatim
from pprint import pprint
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import folium
import networkx as nx
from pyvis.network import Network

class Author:
    """
    Class that represents an author.
    """

    def __init__(self, name, surname):
        self.name = name
        self.surname = surname
        self.collaborators = []
    
    def add_collaborator(self, name, surname):
        self.collaborators.append(Author(name, surname))
    
    def get_collaborators(self):
        ret = []
        for author in self.collaborators:
            ret.append(author.name + ' '  + author.surname)
        
        return ret
    
    def create_df(self):
        pass

class Paper:
    """
    Class that represents a publication.
    """

    def __init__(self, paper):
        self.retrieve_details(paper)
    
    def retrieve_details(self, xml):
        """
        Function that retrieves details of single 
        paper using XML result.
        """

        citation = xml['MedlineCitation']

        self.abstract = str(citation['Article']['Abstract']['AbstractText'][0])

        try:
            self.year = citation['Article']['ArticleDate'][0]['Year']
        except:
            self.year = citation['Article']['Journal']['JournalIssue']['PubDate']['Year']

        self.journal = citation['Article']['Journal']['ISOAbbreviation']
        self.title = citation['Article']['ArticleTitle']

        self.retrieve_authors(xml)
    
    def retrieve_authors(self, xml):
        """
        Function that retrieves all the authors of single paper 
        using XML result.
        """

        citation = xml['MedlineCitation']
        _authors = list(citation['Article']['AuthorList'])
        self.authors = {}

        for _author in _authors:
            affiliations = [aff['Affiliation'] for aff in _author['AffiliationInfo']]
            name = _author['ForeName']
            surname = _author['LastName']
            author = Author(name, surname)
            self.authors[author] = affiliations
    
    def to_dict(self):
        """
        Creates dictionary of information from Paper object.
        """

        return {
            'title': self.title, 
            'year': self.year, 
            'journal': self.journal, 
            'authors': [el.name + ' ' + el.surname for el in self.authors.keys()]
        }

class Search:
    """
    Class that performs search using PubMed with user-defined 
    query.
    """

    def __init__(self, email, retmax, sort):
        self.query = ""
        self.email = email
        self.retmax = retmax
        self.sort = sort
    
    def esearch(self, query):
        """
        Search PubMed based on query.
        """

        self.query = query
        Entrez.email = self.email

        handle = Entrez.esearch(db='pubmed', sort=self.sort, 
                                retmax=self.retmax, retmode='xml', 
                                term=self.query)
        
        return Entrez.read(handle)
    
    def efetch(self, ids_list):
        """
        Fetch titles from PubMed.
        """

        ids = ','.join(ids_list)
        Entrez.email = self.email

        handle = Entrez.efetch(db='pubmed', retmode='xml', 
                                id=ids)
        
        return Entrez.read(handle)
    
    def perform_search(self, query):
        """
        Search and fetch titles from PubMed.
        """

        res = self.esearch(query)
        ids_list = res['IdList']
        
        return self.efetch(ids_list)
    
    def search(self, query):
        """
        Function to actually perform search.
        """

        res = self.perform_search(query)
        papers = []

        for _paper in res['PubmedArticle']:
            try:
                paper = Paper(_paper)
                papers.append(paper)
            except:
                continue
        
        return papers
    
    def create_df(self, papers):
        """
        Creates Pandas DataFrame with information 
        for each retrieved paper based on 
        user-defined query.
        """

        df = pd.DataFrame.from_records([paper.to_dict() for paper in papers])
        df.sort_values(by=['year'], ascending=False, inplace=True)
        df.reset_index(inplace=True, drop=True)

        return df

class Networks:
    """
    Class that creates networks of collaborators.
    """

    def generate_networks_from_papers(self, papers):
        """
        Generates networks of collaborators starting from retrieved papers.
        """

        _df = pd.DataFrame.from_records([paper.to_dict() for paper in papers])
        _collaborators = {}

        for idx in _df.index:
            for _author in _df.loc[_df.index[idx], 'authors']:
                if _author not in _collaborators.keys():
                    _collaborators[_author] = []
                    _collaborators[_author].extend(_df.loc[_df.index[idx], 'authors'])
                else:
                    _collaborators[_author].extend(_df.loc[_df.index[idx], 'authors'])
        
        for _author in _collaborators.keys():
            _collaborators[_author].remove(_author)
            _collaborators[_author] = set(_collaborators[_author])
        
        # Graph based on collaborators for each author
        g = nx.DiGraph()
        g.add_nodes_from(_collaborators.keys())
        for k, v in _collaborators.items():
            g.add_edges_from(([(k, t) for t in v]))
        
        # Interactive network based on graph
        net = Network(height=700, width=700, notebook=True)
        net.from_nx(g)

        return net

class Map:
    """
    Class that creates maps of collaborators.
    """
    
    def extract_geo_from_papers(self, papers):
        """
        Geo-localization based on affiliations authors.
        """

        geo = []

        for _paper in papers:
            for _author in _paper.authors.items():
                for _affiliation in _author[1]:
                    _geo = GeoText(_affiliation).cities
                    if _geo:
                        geo.append({
                            'name': _author[0].name, 
                            'surname': _author[0].surname, 
                            'affiliation': _affiliation, 
                            'geo': ' '.join(_geo)
                        })
        
        return geo
    
    def generate_interactive_map_from_papers(self, geo):
        folium_map = folium.Map()
        geolocator = Nominatim(timeout=10)
        self.num_affiliations_tot = 0
        self.not_found = 0

        for _geo in geo:
            city = _geo['geo']
            self.num_affiliations_tot += 1
            
            try:
                loc = geolocator.geocode(city)
                
                lat = loc.latitude
                lon = loc.longitude
                aff = _geo['affiliation']
                name = _geo['name']
                surname = _geo['surname']
                txt = f'Affiliation: {aff}<br> Name: {name}<br> Surname: {surname}'

                folium.CircleMarker(location=(lat, lon), 
                                    popup=txt, 
                                    fill=True).add_to(folium_map)
            
            except:
                self.not_found += 1
                continue
        
        print(f'Number of total affiliations: {self.num_affiliations_tot}.')
        print(f'Number of affiliations not found: {self.not_found}.')

        return folium_map
