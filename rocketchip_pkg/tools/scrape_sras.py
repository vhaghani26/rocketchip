from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
import time
import sys
import argparse


parser = argparse.ArgumentParser(description='Webscrape NCBI for SRA entries')
parser.add_argument('url', type=str, metavar='<path>', help='Link to url')
arg = parser.parse_args()

driver = webdriver.Firefox()
driver.get(arg.url)


selectables = driver.find_elements(By.NAME, 'EntrezSystem2.PEntrez.Sra.Sra_ResultsPanel.Sra_DisplayBar.Display')

select_for_size = None
for s in selectables:
	if s.text == '20 per page':
		select_for_size = s
		break

select_for_size.click()

time.sleep(1)

max_per_page = driver.find_element(By.ID, 'ps200')
max_per_page.click()

time.sleep(1)

max_page = driver.find_element(By.CLASS_NAME, 'page')
max_page_num = int(max_page.text.split(' ')[2])


for i in range(max_page_num):
	time.sleep(5)
	main_content = driver.find_element(By.ID, 'maincontent')
	reports = main_content.find_elements(By.CLASS_NAME, 'rprt')
	for report in reports: print(report.text)
	next_page = main_content.find_element(By.LINK_TEXT, 'Next >')
	try:
		next_page.click()
	except:
		time.sleep(2)
		popup = WebDriverWait(driver, 10).until(
		EC.presence_of_element_located(
		(By.CLASS_NAME, 'QSIWebResponsiveDialog-Layout1-SI_blp2VywNohCsrQy_button-container')))
		no_thanks = popup.find_elements(By.TAG_NAME, 'button')[1]
		no_thanks.click()
		time.sleep(2)
		next_page.click()

driver.quit()

