	# -*- coding: utf-8 -*-
"""
Created on Mon Oct 11 10:02:37 2021

@author: chrissierwin
"""

import numpy as np

n = 27
user = np.random.randint(0,n)
# add loop with user input buttons that say yes and no
age = str(input('Are you over 21? Y/N '))
age = age.lower() # makes response lower case
if age == 'y' or age == 'yes':
	over_21 = True
else:
	over_21 = False

# Movie and meal
if user == 0:
    movie = 'The Princess and the frog'    # New Orleans, Loiusiana, United States
    meal = 'Red beans & rice or Jambalaya'
    side = 'Cornbread or Cajun Corn Maque Choux '
    drink = 'Pomegranate juice or Sazerac'
    adult_drink = 'Sazerac'
    dessert = "Beignets, King Cake, or Banana's Foster"
    origin = 'American'
    year = '2009'
    
elif user == 1: 
    movie ='Luca'        								# Italy
    meal = 'Lasagna or Linguine pasta'
    side = 'Salad or Garlic bread'
    drink = 'Lemonade'
    adult_drink = 'Limoncello'
    dessert = 'Gelato'
    origin = 'Italian'
    year = '2021'
    
elif user == 2:
    movie = 'Aladdin'									# Middle East, Arabian
    meal = 'Kabobs or Schwarma (chicken)'
    side = 'Basmati rice and sauce'
    drink = 'Limonana'
    adult_drink = 'Arak'
    dessert = 'Monkey bread'
    origin = 'Arabian'
    year = '1992'
        
elif user == 3:
    movie = 'Mulan'										# Japan
    meal = 'Black pepper chicken'
    side = 'Noodles or white rice'
    drink = 'Tea'
    drink_adult = 'Sake'
    dessert = 'Miniature cakes'
    origin = 'Chinese'
    year = '1998'
    
elif user == 4:
    movie = 'Atlantis'
    meal = 'Sheppards pie'
    side = 'Baked beans with bacon or Pasta salad'
    drink = 'Ocean water'
    adult_drink = 'Ocean water jello shots'
    dessert = 'Blue rock candy'
    origin = 'Atlantic Ocean'
    year = '2001'
    
elif user == 5:
    movie = 'The Lion King'								# Africa
    meal = 'Barbecue'
    side = 'Green beans with bacon and Baked beans'
    drink = 'Sparkling fruit punch'
    drink_adult = 'Hakuna Matata cocktail'
    dessert = 'Worms in dirt (ice cream or brownie)'  
    origin = 'African' # based on kenya and tanzania
    year = '1994'
    
elif user == 6:
    movie = 'The Aristocats'							# Paris, France
    meal = 'Steak frites, Croque monsieur, or Boeuf bourgiugnon'
    side = 'Berry spinich salad or French onion soup'
    drink = 'Berry tea'
    adult_drink = 'Cognac or Champagne'
    dessert = 'Cream puffs or Berry tart'  
    origin = 'French'
    year = '1970'
    
elif user == 7:
    movie = 'Raya and the Last Dragon'  # viatnamese
    meal = 'Bahn mi sandwich'
    side = 'Morning glory stir fry or steamed vegetables'
    drink = 'Nuoc mia (sugarcane juice)'
    adult_drink ='Bia hoi (draught beer)'
    dessert = 'Mooncake or Chè'
    origin = 'Asian' # fantasy land of Kumandra, based on southeast asian cultures
    year = '2021'
    
elif user == 8:
    movie = 'A Goofy Movie'					# United States
    meal = 'Meatball sub or Alphabet soup'
    side = 'Potato chips and pickle spears'
    drink = 'Fruit Punch'
    adul_drink ='Spiked punch bowl'
    dessert = 'Red velvet cupcakes'
    origin = 'American'
    year = '1995'
    
elif user == 9:
    movie = 'Coco'									# Mexico
    meal = 'Quesadilla or Tacos'
    side = 'Caprese salad or chips & salsa'
    drink = 'Frozen strawberry lemonade'
    adult_drink = 'Strawberry frozen margarita'
    dessert = 'Tres leches'
    origin = 'Mexican'
    year = '2017'
    
elif user == 10:
    movie = 'Moana'								# Hawaii, United States
    meal = 'Pork hawaiian sliders'
    side = 'Fresh fruit'
    drink = 'Coconut water'
    adult_drink = 'Fish bowl'
    dessert = 'Sweet bread'
    origin = 'Polynesian'
    year = '2016'
    
elif user == 11:
    movie = 'Brave'
    meal = 'Roast beef'
    side = 'Mashed potatoes with gravy and peas'
    drink = 'Spiced cider'
    adult_drink = 'Mead'
    dessert = 'Pie'
    origin = 'Scottish'
    year = '2012'
    
elif user == 12:
    movie = 'Winnie the Pooh'						# London, England
    meal = 'BBQ Chicken Tenders'
    side = 'Honey glazed carrots'
    drink = 'Honey sweet tea'
    adult_drink = "The Bee's Knees cocktail" # https://cookieandkate.com/bees-knees-cocktail-recipe/
    dessert = 'Chocolate pie'
    origin = 'English' # london
    year = '2011'
    
elif user == 13:
    movie = 'Pocahontas'							# United States
    meal = 'Buffalo chicken tenders'
    side = 'Roasted corn and biscuits'
    drink = 'Grape juice'
    adult_drink = 'Red wine'
    dessert = 'Pumpkin pie'
    origin = 'American' # tsenacommacah, virginia
    year = '1995'
    
elif user == 14:
    movie = 'Lady and the Tramp'
    meal = 'Spaghetti and meatballs'
    side = 'Breadsticks'
    drink = 'Italian soda or coffee'
    adult_drink = 'Champagne'
    dessert = 'Tiramisu'
    origin = 'American' # missouri
    year = '1955'
    
elif user == 15:
	movie = 'Hercules'
	meal = 'Gryos'
	takeout = 'Greek takeout'
	side = 'Greek salad and olives'
	drink = 'Frappé'
	adult_drink = 'Retsina'
	dessert = 'Baklava or Lokma'
	origin = 'Greek'
	year = '1997'
	
elif user == 16:
	movie = 'Finding Nemo'
	meal = 'Sushi or Potstickers'
	side = 'Edamame'
	drink = 'Fizzy blue gatorade'
	adult_drink = 'Fizzy blue gatorade with malibu'
	dessert = 'Ice cream filled taiyaki'
	origin = 'Australian'
	year = '2003'
	
elif user == 17:
    movie = 'Tangled'
    meal = 'Pasta with pesto'
    side = 'Spring salad'
    drink = 'Sparkling raspberry lemonade'
    adult_drink = 'Blackberry gin'
    dessert = 'Scones or Blackberry tart'
    origin = 'Normandic'
    year = '2010'
	
elif user == 18:
    movie = 'Chicken Little'
    meal = 'Chicken nuggets'
    side = 'Mashed potatoes and corn'
    drink = 'Orange soda'
    adult_drink = 'Spike orange dreamsicle float'
    dessert = 'Orange dreamsicle float'
    origin = 'American'
    year = '2005'
	
elif user == 19:
    movie = 'Ratatouille'
    meal = 'Ratatouille'
    side = 'Charcuterie or soup'
    drink = 'Grape juice'
    adult_drink = 'Wine'
    dessert = 'Berry tart'
    origin = 'French'
    year = '2007'
	
elif user == 20:
    movie = 'Lilo and Stitch'
    meal = 'Kalua pork sandwich'
    takeout = 'Peanut butter and jelly sandwich with carrot sticks'
    side = 'Roasted pineapple or Roasted carrots'
    drink = 'Ocean water'
    adult_drink = 'Mai tai'
    dessert = 'Açai bowl'
    origin = 'American'
    year = '2002'
			
elif user == 21:
    movie = 'Cars'
    meal = 'Burgers'
    side = 'French fries or Tater tots'
    drink = 'Soda'
    adult_drink = 'Ka-chow cocktail'
    dessert = 'Muddy tracks ice cream'
    origin = 'American'
    year = '2006'
	
elif user == 22:
    movie = 'Pirates of fhe Caribbean'
    meal = 'Jerk chicken'
    side = 'Coconut rice and plantains'
    drink = 'Coconut water or Ocean water'
    adult_drink = 'Jack and coke'
    dessert = 'Rum cake'
    origin = 'Caribbean'
    year = '2003'
	
elif user == 23:
    movie = 'Onward'
    meal = 'Beef stew with rosemary'
    side = 'Steamed green beans'
    drink = 'Juice'
    adult_drink = 'Blackberry rum with edible glitter'
    dessert = 'Strawberry tart with edible glitter'
    origin = 'American'
    year = '2020'
	
elif user == 24:
    movie = 'Incredibles'
    meal = 'Meatloaf'
    side = 'Mashed potatoes and peas'
    drink = 'Coca-cola'
    adult_drink = 'Jack and coke'
    dessert = 'Cake'
    origin = 'American'
    year = '2004'

elif user == 25:
    movie = 'Incredibles 2' 
    meal = 'Chinese take-out'
    side = 'Spring rolls'
    drink = 'Milk tea'
    adult_drink = 'Sake'
    dessert = 'Chinese sesame cookies'
    origin = 'American'
    year = '2018'

elif user == 26:
    movie = 'Toy Story'
    meal = 'Pizza'
    side = 'French fries'
    drink = 'Soda'
    adult_drink = 'Crown and seltzer'
    dessert = 'Ice cream sundae'
    origin = 'American'
    year = '1995'
    
elif user == 27:
    movie = 'The Jungle Book'
    meal = 'Pork and fried plantain skewers'
    side = 'Black bean and white rice'
    drink = 'Coconut water'
    adult_drink = 'Jungle juice'
    dessert = 'Fruit salad'
    origin = 'Indian'
    year = '1967'

elif user == 28:
    movie = 'The Three Musketeers'
    meal = 'Sandwiches'
    side = 'Potato salad'
    drink = 'Apple juice'
    adult_drink = 'Long island tea'
    dessert = 'Pudding'
    origin = 'French'
    year = '2004'

else:
    print('try again')
    

print('Movie: {}'.format(movie))
print('Meal: {}'.format(meal))
print('Side: {}'.format(side))
if over_21 == True:
	print('Alcoholic beverage: {}'.format(adult_drink))
else:
	print('Drink: {}'.format(drink))
print('Dessert: {}'.format(dessert))

'''
# Template
elif user == 2:
    movie = ''
    meal = ''
    side = ''
    drink = ''
    dessert = ''
'''