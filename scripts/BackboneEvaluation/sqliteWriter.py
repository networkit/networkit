import sqlite3
import os

class SqliteResultWriter():

	def __init__(self, dbFile, properties):
		self.dbFile = dbFile

		#Create empty databse
		if not os.path.isfile(dbFile):
			db = sqlite3.connect(dbFile)
			self.createTables(db, properties)
			db.close()

	def createTables(self, db, properties):
		#Graphs table
		db.execute('''CREATE TABLE graphs (name text primary key)''')

		#Algorithms table
		db.execute('''CREATE TABLE algorithms (name text primary key)''')

		#Properties table
		db.execute("CREATE TABLE properties (name text primary key)")

		#Properties table
		db.execute("CREATE TABLE data (graph text, algorithm text, targetEdgeRatio real, property text, value text)")

		db.commit()

	def createRowIfNeccessary(self, db, table, name):
		cursor = db.cursor()
		cursor.execute('''SELECT name FROM ''' + table + ''' WHERE name=?''', (name,))
		data = cursor.fetchone()

		if data is None:
			cursor.execute('''INSERT INTO ''' + table + ''' (name) VALUES (?)''', (name,))

	def receiveResult(self, taskResult):
		db = sqlite3.connect(self.dbFile)
		deleted = [] #List of pairs

		for row in taskResult.data:
			#Create entries in graph and algorithm tables
			algorithmId = row['algorithm']
			graphId = row['graph']
			self.createRowIfNeccessary(db, 'algorithms', algorithmId)
			self.createRowIfNeccessary(db, 'graphs', graphId)

			propertyNames = list(set(row.keys()) - {'graph', 'algorithm', 'targetEdgeRatio'})

			for propertyId in propertyNames:
				ratioId = row['targetEdgeRatio']
				self.createRowIfNeccessary(db, 'properties', propertyId)
				db.execute('''DELETE FROM data WHERE graph=? AND algorithm=? AND targetEdgeRatio=? AND property=?''', (graphId, algorithmId, ratioId, propertyId))
				query = "INSERT INTO data (graph, algorithm, property, targetEdgeRatio, value) VALUES ('" + graphId + "', '" + algorithmId + "', '" + propertyId + "', '" + str(ratioId) + "', '" + str(row[propertyId]) + "')"
				db.execute(query)

		db.commit()
		db.close()
