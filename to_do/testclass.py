class Car:

    def __init__(self, color, model, year):
        self.color = color
        self.model = model
        self.year = year

    def how_old(self, current_year):
        age = current_year - self.year
        return(age)

first_car = Car("blue", "Ranger", 1979)
second_car = Car("red", "Camry", 1992)

first_car.model
second_car.year

first_car.how_old(2023)
second_car.how_old(2023)