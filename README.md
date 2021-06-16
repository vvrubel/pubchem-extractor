PubChem has a standard time limit of 30 seconds per web service request. If a request is not completed within the 30-second limit for any reason, a timeout error will be returned.  To work around certain slower operations, one may use an ‘asynchronous’ approach, where a so-called ‘key’ is returned as a response to the initial request.  This key is then used to check periodically whether the operation has finished, and, when complete, retrieve the results.

All PubChem web pages (or requests to NCBI in general) have a policy that users should throttle their web page requests, which includes web-based programmatic services.  Violation of usage policies may result in the user being temporarily blocked from accessing PubChem (or NCBI) resources.  The current request volume limits are:

- No more than 5 requests per second.
- No more than 400 requests per minute.
- No longer than 300 second running time per minute.

It should be noted that these limits can be lowered through the dynamic traffic control at times of excessive load.  Throttling information is provided in the HTTP header response, indicating the system-load state and the per-user limits.  Based on this throttling information, the user should moderate the speed at which requests are sent to PubChem.

To help maximize uptime and request handling speed, PubChem web servers employ a dynamic, web-request throttling approach that enforces usage policies.  Importantly, during periods of excessive demand, these policies may be dynamically changed to maintain accessibility to all users.  Requests exceeding limits are rejected (HTTP 503 error).  If the user continuously exceeds the limit, they will be blocked for a period of time.