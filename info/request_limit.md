- PubChem has a standard `time limit of 30 seconds per web service request`. 
  If a request is not completed within the 30-second limit for any reason, a timeout error will be returned.
  
- To work around certain slower operations, one may use an `asynchronous` approach, 
  where a so-called `key` is **returned as a response** to the initial request.
  This `key` is used to **check** periodically whether **the operation has finished**, and, when complete, retrieve the results.

- All PubChem web pages (or requests to NCBI in general) have a policy that users should `throttle their web page requests`, 
  which includes web-based programmatic services.  
  Violation of usage policies may result in the user **being temporarily blocked** from accessing PubChem (or NCBI) resources.  

> The current request volume limits are:
> 
- No more than 5 requests per second.
- No more than 400 requests per minute.
- No longer than 300 second running time per minute.

It should be noted that these limits can be lowered through the dynamic traffic control at times of excessive load.  

`Throttling information` is provided in the `HTTP header response`, indicating the system-load state and the per-user limits.  
Based on this throttling information, the user should moderate the speed at which requests are sent to PubChem.

---

To help maximize uptime and request handling speed, PubChem web servers employ a dynamic, web-request throttling approach that enforces usage policies.  
Importantly, during periods of excessive demand, these policies may be dynamically changed to maintain accessibility to all users.  
Requests exceeding limits are rejected (HTTP 503 error).  
If the user continuously exceeds the limit, they will be blocked for a period of time.


Therefore, the user should moderate the speed at which requests are sent to PubChem, according to the traffic status of PubChem 
and the extent to which the user is approaching limits.  
This information is provided in specialized HTTP response headers accompanying all PUG-REST web requests.  
For example, the HTTP response header contains a line similar to the following:
```
X-Throttling-Control: Request Count status: Green (0%), Request Time status: Green (0%), Service status: Green (20%)
```

The first two status indicators (`Request Count` and `Time statuses`) give information on your usage of the service in one of four states:
```
Green - less than 50% of the permitted request limit has been used
Yellow - between 50% and 75% of the request limit has been used
Red - more than 75% of the request limit has been reached
Black - the limit has been exceeded and requests are being blocked
```
The third indicator (Service status) shows the concurrent usage of the service in one of four states:
```
Idle (Green) - Low concurrent usage being applied to the service at present
Moderate (Yellow) - a moderate number of concurrent requests are being handled
Busy (Red) - a significant number of concurrent requests are being handled
Overloaded (Black) - an excessively high number of concurrent requests are being handled
```
---

It is important to note that there are many instances of PubChem services running in parallel. 
Each instance receives traffic from a load balancer, which distributes the requests across the system. 
Thus, when a stream of requests is sent to PubChem, the responses will be relative to the PubChem server instance handling the request. 

One server instance can become overloaded while others may not, depending on the overall nature of requests sent to that server. 
When providing many requests, one should moderate the speed requests are sent to according to the worst-case usage feedback received.  
This will prevent uneven rejection of requests by PubChem services.

- Downloads through the `structure download tool` are limited to a `maximum of 500 000 records per request`, 
  with an exception of `image downloads`, which are `limited to 50 000 per request`. 
  These limits keep the download file sizes practical. 
  
- `The structure download service` is accessible using the `PubChem Power User Gateway (PUG)`, which allows for a programmatic access to PubChem. 
  One can download more than this limit, through multiple interactive or programmatic requests.