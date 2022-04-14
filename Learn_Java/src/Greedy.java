// FROM: https://www.geeksforgeeks.org/activity-selection-problem-greedy-algo-1/
// the following implementation assumes that the activities
// are already sorted according to their finish time
import java.util.*;
import java.lang.*;
import java.io.*;

public class Greedy{
    // prints a max set of activities that can be done by a 
    // single person one at a time
    // n:   total number of activities
    // s[]: an array that contains start time of all activities
    // f[]: an array that contains finish time of all activities
    public static void printMaxActivities(int s[], int f[], int n){
        int i, j;
        System.out.println("The following activities are selected: n");

        // the first activity gets selected
        i = 0;
        System.out.println(i + " ");

        // consider the rest of the activities
        for (j=1; j<n; j++){
            // if this activity has a start time <= finish time of the 
            // previously selected activity, then select it
            if (s[j] >= f[i]){
                System.out.println(j + " ");
                i = j;
            }
        }
    }
    public static void main(String[] args){
        int s[] = {1, 3, 0, 5, 8, 5};
        int f[] = {2, 4, 6, 7, 9, 9};
        int n = s.length;

        printMaxActivities(s, f, n);
    }
}

// public class Greedy{
//     private boolean isCounterEnabled = true;
//     private int counter = 4;

//     @Getter @ Setter
//     List users;

//     public Greedy(){
//         users = new ArrayList<>();
//     }

//     public boolean switchCounter(){
//         this.isCounterEnabled = !this.isCounterEnabled;
//         return this.isCounterEnabled;
//     }

// }

