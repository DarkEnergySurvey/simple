import urllib
import yaml
import os
from PIL import Image

with open('config.yaml', 'r') as ymlfile:
    cfg = yaml.load(ymlfile)

save_dir = os.path.join(os.getcwd(), cfg['output']['save_dir'])
if not os.path.exists(save_dir):
    os.mkdir(save_dir)

def retrieve_des_dr1_image(filename, ra, dec):
    filename = '{}/image_{:0.2f}_{:0.2f}.png'.format(save_dir, ra, dec)
    url = "http://legacysurvey.org/viewer/jpeg-cutout?ra={0}&dec={1}&zoom={2}&layer=des-dr1".format(ra, dec, 10)
    urllib.urlretrieve(url, filename) #Retreaves and saves each image

    return(url)

def retrieve_ps1_image(filename, ra, dec):
    filename = '{}/image_{:0.2f}_{:0.2f}.png'.format(save_dir, ra, dec)
    #page_url = "https://ps1images.stsci.edu/cgi-bin/ps1cutouts?pos={0}%2C%2B{1}&filter=g&filetypes=stack&auxiliary=data&size={2}&output_size=512&verbose=0&autoscale=90.0&catlist=".format(ra, dec, 6000)
    page_url = "https://ps1images.stsci.edu/cgi-bin/ps1cutouts?pos={0}%2C{1}&filter=g&filetypes=stack&auxiliary=data&size={2}&output_size=512&verbose=0&autoscale=90.0&catlist=".format(ra, dec, 6000)
    image_line = [line for line in urllib.urlopen(page_url).read().split("\n") if "<td><img src=" in line][0]
    url = "http://" + image_line[image_line.index("/")+2:].split("\"", 1)[0]
    urllib.urlretrieve(url, filename) #Retreaves and saves each image

    return(page_url)

def retrieve_image(filename, ra, dec, survey):
    if(survey == 'des'):
        return(retrieve_des_dr1_image(filename, ra, dec))
    elif(survey == 'ps1'):
        return(retrieve_ps1_image(filename, ra, dec))
    else:
        return("No image found")
