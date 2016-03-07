# @Author 
# Christine Staiger
# staiger@cwi.nl; staigerchristine@gmail.com
# July 2013
from couchdb import Server
server = Server('https://picas.grid.sara.nl:6984')
db = server['aces']
for x in db:
    print x

#deleting items forever:
db.__delitem__(x['_id'])

from picas.clients import CouchClient
client = CouchClient(url="https://picas.grid.sara.nl:6984", db="aces", username="cstaiger", password="dq4^#sdk")
db = client.db
for x in db:
    doc = db[x]
    db.(doc)

docs = []
for x in db:
    doc = db[x]
    docs.append(doc)

for d in docs:
    d['lock'] > 0
    print d['output']


