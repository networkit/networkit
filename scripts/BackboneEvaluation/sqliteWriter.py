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
		query = '''CREATE TABLE properties ('''
		for p in properties:
			for key in list(p.getTypes().keys()):
				query += key + ' ' + p.getTypes()[key] + ', '
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
		#db.execute('''UPDATE graphs SET name=?, loadingTime=? WHERE id=?''', (graphName, taskResult.loadingTime, graphId))

		for row in taskResult.data:
			#Create entries in graph and algorithm tables
			algorithmId = row['algorithm']
			graphId = row['graph']
			self.createRowIfNeccessary(db, 'algorithms', algorithmId)
			self.createRowIfNeccessary(db, 'graphs', graphId)

			#Delete existing properties
			db.execute('''DELETE FROM properties WHERE graph=? AND algorithm=?''', (graphId, algorithmId))

			#Insert datarow into database
			propertyNames = list(row.keys())
			query = "INSERT INTO properties ("
			for propertyName in propertyNames:
				query += propertyName + ', '
			query = query[:-2] + ") VALUES ("
			for propertyName in propertyNames:
				query += "'" + str(row[propertyName]) + "', "
			query = query[:-2] +");"
			db.execute(query)

		db.commit()
		db.close()
