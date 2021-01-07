#!/usr/bin/env python
# coding: utf-8

# In[1]:


import math as m
import pandas as pd

#declare global constants
pi = m.pi

a = 6378137
f_inv = 298.257223563

#flattening
f = pow(f_inv,-1)

#(4.9)
b = a - f*a

#first eccentricity e^2 (4.10)
e_2 = ( a**2 - b**2 ) / a**2

#second eccentricity e'2 (4.11)
e__2 = ( a**2 - b**2 ) / b**2


class Point():
    """
    Description:
        This class is set up to handle one set of coordinated being 
        input at a time and to convert them appropriatly
    """
    
    def __init__(self, one, two, three, entered = "Curvlinear", four=-1):
        """
        Description:
            Initializes all needed values within the point. Takes one set of 
            values and initializes the class variables of the other type of 
            coordinates
            
        Input:
            one, is either O or x
            two, is either L or y
            three, is either h or z
            entered = "Curvlinear" or "Cartesian"
            
        Output:
        """
        #counter for stuff
        self.count = 0
        
        if entered == "Curvlinear":
            #declaration style one
            self.O = m.radians(one)
            self.L = m.radians(two)
            #H + N
            self.h = three + four
            
            self.curv_to_cart()
            
        elif entered == "Cartesian":
            #declaration style two
            self.X = two #longitude
            self.Y = one #latitude
            self.Z = three
            
            self.cart_to_curv()
            
        else:
            print("Inproper declaration")
           # raise
            
    def curv_to_cart(self):
        """
        Description:
            Converts curvlinear coordinates to cartesian coordinates

        Input:
            
        Output:
            self.X, self.Y, and self.Z are declared
        """
        N = a / (m.sqrt(1 - e_2 * m.sin(self.O)**2))
        
        self.X = (N + self.h)* m.cos(self.O) * m.cos(self.L)
        self.Y = (N + self.h) * m.cos(self.O) * m.sin(self.L)
        self.Z = ((1 - e_2) * N + self.h) * m.sin(self.O)

    def cart_to_curv(self):
        """
        Description:
            Converts cartesian coordinates curvlinear coordinates

        Input:

        Output:

        """
        #(4.54 - 4.59)

        #these will remain the same
        self.L = m.atan2(self.Y,self.X)
        self.p = m.sqrt(self.X**2 + self.Y**2)

        #initialize height, latitude, prime vertical radius
        h_0 = 0
        O_0 = m.atan(self.Z / self.p * 1 / (1 - e_2))
        
        convergence = .01

        #begin iterating
        self.cart_to_curv_converge(convergence, O_0, h_0)
        
    
    def cart_to_curv_converge(self, convergence, O_k_1, h_k_1):
        """
        Description:
            Repetitively iterates this function until the convergence parameter 
            is reached

        Input:
            convergence, the decimal to monitor for convergence reached
            O_k_1, the value of the past iteration
            h_k_1, the value of the past iteration

        Output:

        """
        #(4.54 - 4.59)
        #kth vertical radius
        N_k = a / m.sqrt(1 - e_2 * pow(m.sin(O_k_1),2))

        #new height
        h_k = ( self.p / m.cos(O_k_1) ) - N_k
        
        #new latitude
        O_k = m.atan( self.Z / self.p * pow( (1-e_2 * N_k/(N_k + h_k)), -1))
                
        if (abs(h_k_1 - h_k) < convergence):
            #convergence has been met and the values can be declared
            self.O = O_k
            self.h = h_k
            return

        #convergence has not been met, must continue to iterate
        return self.cart_to_curv_converge(convergence, O_k, h_k)
    
    
class Coordinates():
    """
    Contains a list of points with reading and writing functions for them
    """
    
    def __init__(self, file, entered = "Curvlinear"):
        """
        Description:
        
        Input:
            file, 
            entered = "Curvlinear"
            
        Output:
            
        """
        self.file = file
        self.points = []
        
        self.initialization(entered)
        
        #initialize this dataframe
        self.coords_to_df()
        
    def initialization(self, entered):
        """
        Description:
            Initializes the curvlinear variation of reading in
            
        Style: point#   one   two   three
        """
        names = ['one','two','three']
        if entered == "Curvlinear":
            names = ['one','two','three','four']
            
        df = pd.read_csv(self.file, 
                        index_col=0, 
                        sep="\t", 
                        header = None,
                        names = names) 
        
        for index, row in df.iterrows():
            if entered == "Curvlinear":
                self.points.append(Point(row['one'],
                                    row['two'],
                                    row['three'],
                                    "Curvlinear",
                                    row['four']))
            elif entered == "Cartesian":
                self.points.append(Point(row['one'],
                                    row['two'],
                                    row['three'],
                                    "Cartesian"))
            
    def coords_to_df(self):
        """
        Description:
            return a dafaframe of all values
            
            rows "O", "L", "h". "X", "Y", "Z"
        """
        self.df = pd.DataFrame(columns = ["O", "L", "h", "X", "Y", "Z"])
        
        for point in self.points:
            self.df = self.df.append({"O": m.degrees(point.O), 
                                    "L": m.degrees(point.L), 
                                    "h": point.h,
                                    "X": point.X, 
                                    "Y": point.Y, 
                                    "Z": point.Z},  
                                    ignore_index = True) 
            
    def straight_line_distance(self, i_1, i_2):
        """
        Description:
            Finds the distances between the two points
            
        Input:
            The indexes of the points in self.points to find the distance from
            
        Output:
            Returns a value in meters
        """
        Xd = self.points[i_1].X - self.points[i_2].X
        Yd = self.points[i_1].Y - self.points[i_2].Y
        Zd = self.points[i_1].Z - self.points[i_2].Z
        
        return m.sqrt(Xd**2 + Yd**2 + Zd**2)
    
    def convert_back(self):
        """
        Description:
            After the spherical coordinates were put it, calling this command 
            will back convert the coords via the X, Y, Z ones that were calculated
            
        Input:
            none
            
        Output:
            Updates all point O, L, h values
            Updates the coordinates' datafram self.df
        """
        #back converting O, L, h
        #converts each point
        for point in self.points:
            point.cart_to_curv()

        #updates the dataframe within the coords object
        self.coords_to_df()


# In[2]:


#GPSleveling_CANADA_Curv.txt
#GPSleveling_CANADA.txt
coords = Coordinates('GPSleveling_CANADA.txt', entered = "Curvlinear")


# In[3]:


#without back convertion check
coords.df


# In[5]:


#back converting O, L, h
coords.convert_back()

#back converted coords
coords.df


# In[6]:


#calgary and honolulu coords
coords = Coordinates('GPSleveling_CANADA_Curv.txt', entered = "Curvlinear")

#without back convertion check
coords.df


# In[7]:


#back converting O, L, h
coords.convert_back()

#back converted coords
coords.df


# In[8]:


d = coords.straight_line_distance(0,1)
print("Straight line distance is: " + str(d) + "m")


# In[ ]:




