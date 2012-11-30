import os
import sys
import re
import gc
import psutil
import time
from datetime import datetime

import sqlite3
#import MySQLdb

#Reading the arguments
Argument = []
Argument = sys.argv[1:] 
app_folder = Argument[0]


def Scheduler():

	try:
		#Connect to database
		#conn = MySQLdb.connect (host = "localhost", user = "mdat", passwd = "mdat", db = "mdat")
		db_path = os.path.join(app_folder, 'databases', 'storage.db')
		con = sqlite3.connect(db_path)
		con.row_factory = sqlite3.Row
		c=con.cursor()

		#Check lock
		sys_info_sql = 'SELECT * FROM systems WHERE id = "1"'
		c.execute(sys_info_sql)

	except:
		print "Error (1) with database in Scheduler:", sys.exc_info()[0]
		con.close()
		gc.collect()
		return 0
		
	row = c.fetchone()
	locked = row['locked']
	default_process = row['default_process']
	max_process = row['max_process']

	if (locked == 'F'):	#Lock is free: run analysis
		try:
			#Lock
			locked_sql = 'UPDATE systems SET locked = "T" WHERE id = "1"'
			with con:
				c.execute(locked_sql)
		except:
			print "Error (2) with database in Scheduler:", sys.exc_info()[0]
			con.close()
			gc.collect()
			return 0

		#Read datasets
		c.execute('SELECT id, status FROM results ORDER BY id ASC')
		results = c.fetchall()

		#Count the number of running processes
		running_process=0
		for result in results:
			if (result['status']=='Running'):
				running_process = running_process + 1

		#Running process < default, continue
		if (running_process < default_process):

			remaining_process=0

			for result in results:

				if (result['status']=='Waiting'):

					#Take dataset id
					result_id = result['id']

					#Start, set status to Running
					try:
						start_time = datetime.now()
						
						result_update_sql = 'UPDATE results SET status = "Running", start_time = "' + str(start_time) + '" WHERE id = "' + str(result_id) + '"'
						with con:
							c.execute(result_update_sql)

					except:
						print "Error (3) with database in Scheduler:", sys.exc_info()[0]
						con.close()
						gc.collect()
						break

					remaining_process = remaining_process + 1

					try:
						executer_path = os.path.join(app_folder, 'modules', 'Executer.py')
						cmd = 'python ' + executer_path + ' ' + app_folder + ' ' + str(result_id) + ' &'
						os.system(cmd)
					except:
						print "Error (3) with analyzing process in Scheduler:", sys.exc_info()[0]
						con.close()
						gc.collect()
						break

					#Check number of remaining allowed processes
					if (remaining_process >= default_process - running_process):
						avail_phymem = psutil.avail_phymem()
						avail_virtmem = psutil.avail_virtmem()
						if (avail_phymem < 500000000) | (avail_virtmem < 1000000000) | (remaining_process >= max_process - running_process):
							break

		#Unlock and Halt
		try:
			unlocked_sql = 'UPDATE systems SET locked = "F" WHERE id = "1"'
			with con:
				c.execute(unlocked_sql)

		except:
			print "Error (4) with database in Scheduler:", sys.exc_info()[0]
			con.close()
			gc.collect()
			return 0

		con.close()
		gc.collect()
		return 0

	else:	#Lock is not free: halt
		con.close()
		gc.collect()
		return 0

#Main program
if __name__ == "__main__":
	Scheduler()
