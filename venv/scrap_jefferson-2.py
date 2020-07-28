from selenium import webdriver
import time
from bs4 import BeautifulSoup
import lxml
from webdriver_manager.chrome import ChromeDriverManager
global driver
global output
def main():
    get_off_targets()

def Actual_genomic_hit(output):

    rows = BeautifulSoup(driver.page_source, 'lxml').find(id='myTable0').tbody.find_all('tr')[0:5]
    for i in range(5):
        targets =[]
        targets.append(rows[i].find("td", {"style": "white-space: nowrap;font-family:Courier New;"}).text)
    output[i] = list


def get_off_targets():
    options = webdriver.ChromeOptions()
    options.add_argument('headless')
    options.add_argument('window-size=1920x1080')
    options.add_argument("disable-gpu")
    global driver
    driver = webdriver.Chrome(ChromeDriverManager().install(), options=options)
    driver.set_page_load_timeout(20)


    output = {}
    with open('just_guides.txt', 'r') as inputs_file:
        lines = inputs_file.readlines()
        for line in lines:
            if line != '\n' and line != '':
                driver.get('https://cm.jefferson.edu/Off-Spotter/')
                time.sleep(5)
                print(line)
                driver.find_element_by_id('input20mers').send_keys(line)
                driver.find_element_by_id('submitmerSearch').click()
                time.sleep(5)
                Actual_genomic_hit(output)
    driver.close()
    print(output)
main()