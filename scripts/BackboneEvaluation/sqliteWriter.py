import sqlite3
import os

class SqliteResultWriter():

	def __init__(self, dbFile):
		self.dbFile = dbFile

		#Create empty databse
		if not os.path.isfile(dbFile):
			db = sqlite3.connect(dbFile)
			self.createTables(db)
			db.close()

	def createTables(self, db):
		#Graphs
		db.execute('''CREATE TABLE graphs (
			id integer primary key autoincrement,
			name text,
			loadingTime real)''')

		#Algorithms
		db.execute('''CREATE TABLE algorithms (
			id integer primary key autoincrement,
			name text)''')

		#Properties
		db.execute('''CREATE TABLE properties (
			graphId integer,
			algorithmId integer,
			clusteringCoefficient real,
			cpvDistanceFromOriginal real,
			cpvDistanceFromOriginalNormalized real,
			degreeDistCoefficient real,
			diameter integer,
			graphStructuralRandMeasure real,
			keptEdgesPercent real,
			keptNodesPercent real,
			numEdges integer,
			numNodes integer,
			runningTime real)''')

	def createRowIfNeccessary(self, db, table, name):
		cursor = db.cursor()
		cursor.execute('''SELECT id FROM ''' + table + ''' WHERE name=?''', (name,))
		data = cursor.fetchone()

		if data is None:
			cursor.execute('''INSERT INTO ''' + table + ''' (name) VALUES (?)''', (name,))
			return cursor.lastrowid
		else:
			return data[0]


	def receiveResult(self, taskResult):
		db = sqlite3.connect(self.dbFile)

		#Graph
		graphName = taskResult.task.graphPath
		graphId = self.createRowIfNeccessary(db, 'graphs', graphName)
		db.execute('''UPDATE graphs SET name=?, loadingTime=? WHERE id=?''', (graphName, taskResult.loadingTime, graphId))

		for backbone in taskResult.backboneProperties:
			#Algorithms
			algorithmName = backbone.algorithmName
			algorithmId = self.createRowIfNeccessary(db, 'algorithms', algorithmName)
			db.execute('''UPDATE algorithms SET name=? WHERE id=?''', (algorithmName, algorithmId))

			#Delete existing properties
			db.execute('''DELETE FROM properties WHERE graphId=? AND algorithmId=?''', (graphId, algorithmId))

			#Properties
			propertyNames = ["clusteringCoefficient", "cpvDistanceFromOriginal",
				"cpvDistanceFromOriginalNormalized", "degreeDistCoefficient", "diameter",
				"graphStructuralRandMeasure", "keptEdgesPercent", "keptNodesPercent",
				"numEdges", "numNodes", "runningTime"]
			values = [graphId, algorithmId] + list(map(lambda x: getattr(backbone, x), propertyNames))
			db.execute('''INSERT INTO properties VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)''', tuple(values))

		db.commit()
		db.close()
