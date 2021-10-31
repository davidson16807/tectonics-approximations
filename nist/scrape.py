import time
import re
import requests

compounds = {
	'Water':              'C7732185',
	'Nitrogen':           'C7727379',
	'Hydrogen':           'C1333740',
	'Oxygen':             'C7782447',
	'Carbon monoxide':    'C630080',
	'Carbon dioxide':     'C124389',
	'Methane':            'C74828',
	'Ethane':             'C74840',
	'Helium':             'C7440597',
	'Argon':              'C7440371',
	'Ammonia':            'C7664417',
	'Tetrafluoromethane': 'C75730',
	'Benzene':            'C71432',
	'Sulfur dioxide':     'C7446095',
}
pressures = ['0.0', '0.03', '0.1', '0.3', '1.0', '3.0', '10.0']

concatenation = ''
for name, id in compounds.items():
	for p in pressures:
		print(f'{name}, {p}MPa')
		text = requests.get(
			f'''https://webbook.nist.gov/cgi/fluid.cgi?Action=Data&Wide=on&ID={id}&Type=IsoBar&Digits=5&P={p}&THigh=&TLow=&TInc=50.0&RefState=DEF&TUnit=K&PUnit=MPa&DUnit=kg%2Fm3&HUnit=kJ%2Fkg&WUnit=m%2Fs&VisUnit=uPa*s&STUnit=N%2Fm'''
		).text.replace('\n', f'\t{name}\n')
		if len(concatenation) > 0:
			text = text.split('\n',1)[1]
		concatenation += text
		time.sleep(5)

with open('test.csv', 'w') as f: f.write(concatenation)
