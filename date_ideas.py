# date ideas

# import packages
from datetime import date
import random

# grab today's date from the computer
month = int(date.today().strftime("%m"))

# date ideas as lists
go_out = ['play putt putt', 'local zoo', 'drive-in movie', 'bowling', 'visit an arcade', 'try an escape room',
          'visit a museum', 'visit an aquarium', 'first friday', 'karaoke night', 'rock climbing', 
          'go hiking', 'visit a new city', 'go to a concert', 'take a dance lesson', 'try a new restaurant',
          'wine tasting tour', 'ice cream date', 'play lasertag', 'play trivia', 'thrift store', "farmer's market",
          'visit a planetarium', 'visit a theme park', 'go to an art gallery', 'go to the beach', 'go to a bar', 
          "get a couple's massage", 'take a mini road trip', 'visit a bookstore', 'visit a local park']

at_home = ['cook something new', 'movie marathon', 'play board games', 'video game tournament', 'make a pillowfort',
           'have a spa night', 'plan your next vacation', 'have a candlelit dinner']

winter = ['sledding', 'build a snowman', 'gingerbread house contest', 'ice skating', 'bake cookies',
          'cuddle up with a movie']

spring = ['picnic', 'play baseball', 'watch the sunset', 'play frisbee', 'go paddleboarding', 'paint at a local park',
          'pick flowers at a local park']

summer = ['go swimming', 'stargaze', 'go to a sporting event', 'have a bonfire and roast marshmallows', 
          'go camping', 'go out for ice cream']

fall = ['go to a pumpkin patch', 'carve pumpkins', 'make caramel apples', 'pumpkin spice things up', 'take a hayride',
        'make a pie']

# suggest date ideas based on the current season
if (month == 12) or (month == 1) or (month == 2):
    season = winter
elif (month == 3) or (month == 4) or (month == 5):
    season = spring
elif (month == 6) or (month == 7) or (month == 8):
    season = summer
elif (month == 9) or (month == 10) or (month == 11):
    season = fall

# let user choose the type of date they'd like
user_input = input('Do you want to look at dates by season, going out, or staying home? ')

# allow user input to match code
if user_input == 'season':
    user_input = season
elif user_input == 'going out':
    user_input = go_out
elif user_input == 'staying home':
    user_input = at_home

# choose some date and print result :)
date_num = len(user_input)
num = random.randint(0,date_num)
print(user_input[num])