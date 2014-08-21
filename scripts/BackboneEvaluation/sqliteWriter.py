import sqlite3
import os

class SqliteResultWriter():

	def __init__(self, dbFile. properties):
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
		query = '''CREATE TABLE properties (graph text, algorithm text, '''
		for p in properties:
			for key, typee in p.getTypes().iteritems():
				query += key + ' ' + typee + ', '
		query = query[:-2] + ")"
		db.execute(query)

	def createRowIfNeccessary(self, db, table, name):
		cursor = db.cursor()
		cursor.execute('''SELECT name FROM ''' + table + ''' WHERE name=?''', (name,))
		data = cursor.fetchone()

		if data is None:
			cursor.execute('''INSERT INTO ''' + table + ''' (name) VALUES (?)''', (name,))

	def receiveResult(self, taskResult):
		db = sqlite3.connect(self.dbFile)

		#Graph
		graphId = taskResult.task.graph.name
		self.createRowIfNeccessary(db, 'graphs', graphId)
		#db.execute('''UPDATE graphs SET name=?, loadingTime=? WHERE id=?''', (graphName, taskResult.loadingTime, graphId))

		for row in taskResult.data:
			algorithmId = data['algorithm']
			self.createRowIfNeccessary(db, 'algorithms', algorithmId)

			#Delete existing properties
			db.execute('''DELETE FROM properties WHERE graph=? AND algorithm=?''', (graphId, algorithmId))

			#Insert datarow into database
			propertyNames = list(row.keys())
			query = "INSERT INTO properties ("
			for propertyName in propertyNames:
				query += propertyName + ', '
			query = query[:-2] + ") VALUES ("
			for propertyName in propertyNames:
				query += row[propertyName] + ', '
			query = query[:-2] +");"
			db.execute(query)

		db.commit()
		db.close()
